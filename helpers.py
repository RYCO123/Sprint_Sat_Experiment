import numpy as np
import requests
import time
from datetime import datetime
from constants import R_EARTH, MU, INITIAL_ALTITUDE, INCLINATION_DEG, CACHE_EXPIRATION_SECONDS

# --- API Caching ---
api_cache = {
    "atmosphere": {"time": 0, "data": None, "last_coords": None},
    "magnetic_field": {"time": 0, "data": None, "last_coords": None},
    "space_weather": {"time": None, "data": (150.0, 150.0, 12)}
}

def calculate_circular_orbit_conditions():
    """
    Calculate initial conditions for a circular orbit at specified altitude.
    
    Returns:
        np.array: Initial state vector [x0, y0, z0, vx0, vy0, vz0]
    """
    r_orbit = R_EARTH + INITIAL_ALTITUDE
    v_circular = np.sqrt(MU / r_orbit)
    
    # Initial position (perigee)
    x0 = r_orbit
    y0 = 0.0
    z0 = 0.0
    
    # Initial velocity (circular orbit with inclination)
    inclination_rad = np.radians(INCLINATION_DEG)
    vx0 = 0.0
    vy0 = v_circular * np.cos(inclination_rad)
    vz0 = v_circular * np.sin(inclination_rad)
    
    return np.array([x0, y0, z0, vx0, vy0, vz0])

def get_geodetic_coords(state):
    """
    Convert ECI state to geodetic coordinates.
    
    Args:
        state: State vector [x, y, z, vx, vy, vz]
    
    Returns:
        tuple: (latitude_deg, longitude_deg, altitude_m)
    """
    r = np.linalg.norm(state[:3])
    alt = r - R_EARTH
    lat_rad = np.arcsin(state[2] / r)
    lon_rad = np.arctan2(state[1], state[0])
    return np.degrees(lat_rad), np.degrees(lon_rad), alt

def calculate_altitude(state):
    """
    Calculate altitude from state vector.
    
    Args:
        state: State vector [x, y, z, vx, vy, vz]
    
    Returns:
        float: Altitude in meters
    """
    return np.linalg.norm(state[:3]) - R_EARTH

def get_space_weather_data(date):
    """
    Fetches F10.7 and Ap solar indices from the CelesTrak API.
    
    Args:
        date: datetime object
    
    Returns:
        tuple: (f107_avg, f107_daily, ap_daily)
    """
    global api_cache
    current_day = date.strftime('%Y-%m-%d')
    if current_day == api_cache["space_weather"].get("time"):
        return api_cache["space_weather"]["data"]
    
    try:
        url = "https://celestrak.org/SpaceData/sw19571001.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        lines = response.text.strip().split('\n')
        
        for line in reversed(lines):
            if line.startswith(f"{date.year} {date.month:02d} {date.day:02d}"):
                parts = line.split()
                ap_daily, f107_daily, f107_avg = int(parts[8]), float(parts[12]), float(parts[13])
                result = (f107_avg, f107_daily, ap_daily)
                api_cache["space_weather"] = {"data": result, "time": current_day}
                print(f"Fetched Space Weather for {current_day}: F10.7_avg={f107_avg}, F10.7_daily={f107_daily}, Ap={ap_daily}")
                return result
    except requests.exceptions.RequestException as e:
        print(f"WARN: Could not fetch space weather data: {e}. Using cached/default values.")
    
    return api_cache["space_weather"]["data"]

def get_atmospheric_density_nasa(sim_datetime, lat, lon, alt_m, f107a, f107, ap):
    """
    Gets atmospheric density by calling the official NASA CCMC API.
    
    Args:
        sim_datetime: datetime object
        lat, lon: Latitude and longitude in degrees
        alt_m: Altitude in meters
        f107a, f107, ap: Space weather parameters
    
    Returns:
        float: Atmospheric density in kg/m³
    """
    global api_cache
    coords = (round(lat, 2), round(lon, 2), round(alt_m / 1000, 2))
    
    if (time.time() - api_cache["atmosphere"]["time"] < CACHE_EXPIRATION_SECONDS and
            api_cache["atmosphere"]["last_coords"] == coords):
        return api_cache["atmosphere"]["data"]
    
    alt_km = alt_m / 1000.0
    url = (f"https://ccmc.gsfc.nasa.gov/api/ws/nrlmsise00/get/+"
           f"year/{sim_datetime.year}/month/{sim_datetime.month}/day/{sim_datetime.day}/+"
           f"hour/{sim_datetime.hour}/min/{sim_datetime.minute}/sec/{sim_datetime.second}/+"
           f"lat/{lat}/lon/{lon}/alt/{alt_km}/"
           f"f107a/{f107a}/f107/{f107}/ap/{ap}")
    
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        data = response.json()['data-request']
        density_g_cm3 = float(data[0]['total_mass_density'])
        density_kg_m3 = density_g_cm3 * 1000.0
        api_cache["atmosphere"] = {"time": time.time(), "data": density_kg_m3, "last_coords": coords}
        return density_kg_m3
    except (requests.exceptions.RequestException, ValueError, KeyError, IndexError) as e:
        cached_data = api_cache["atmosphere"]["data"]
        if cached_data is not None:
            return cached_data
        else:
            return 1e-15

def get_magnetic_field_bgs(sim_datetime, lat, lon, alt_m):
    """
    Gets magnetic field vector by calling the BGS IGRF API.
    
    Args:
        sim_datetime: datetime object
        lat, lon: Latitude and longitude in degrees
        alt_m: Altitude in meters
    
    Returns:
        dict: Magnetic field components in nT {'x': north, 'y': east, 'z': down}
    """
    global api_cache
    coords = (round(lat, 2), round(lon, 2), round(alt_m / 1000, 2))
    
    if (time.time() - api_cache["magnetic_field"]["time"] < CACHE_EXPIRATION_SECONDS and
            api_cache["magnetic_field"]["last_coords"] == coords):
        return api_cache["magnetic_field"]["data"]
    
    alt_km = alt_m / 1000.0
    date_str = sim_datetime.strftime('%Y-%m-%d')
    url = (f"https://geomag.bgs.ac.uk/ws/igrf/?"
           f"latitude={lat}&longitude={lon}&altitude={alt_km}&date={date_str}&format=json")
    
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        data = response.json()
        components = data['geomagnetic-field-model-result'][0]['field-components']
        b_field_ned_nT = {"x": components['north'], "y": components['east'], "z": components['down']}
        api_cache["magnetic_field"] = {"time": time.time(), "data": b_field_ned_nT, "last_coords": coords}
        return b_field_ned_nT
    except (requests.exceptions.RequestException, ValueError, KeyError, IndexError) as e:
        cached_data = api_cache["magnetic_field"]["data"]
        if cached_data is not None:
            return cached_data
        else:
            return {"x": 0, "y": 0, "z": 0}

def calculate_sun_vector(sim_time):
    """
    Calculate sun direction vector (simplified model).
    
    Args:
        sim_time: Simulation time in seconds
    
    Returns:
        np.array: Sun direction vector [x, y, z]
    """
    # Simplified sun position (assume sun is in the x-y plane)
    sun_angle = (sim_time / (365.25 * 86400.0)) * 2 * np.pi
    return np.array([np.cos(sun_angle), np.sin(sun_angle), 0])

def check_eclipse(state, sun_vector):
    """
    Check if satellite is in Earth's shadow.
    
    Args:
        state: State vector [x, y, z, vx, vy, vz]
        sun_vector: Sun direction vector
    
    Returns:
        bool: True if in eclipse, False otherwise
    """
    r_sat = state[:3]
    sun_sat_dot = np.dot(r_sat, sun_vector)
    r_sat_mag = np.linalg.norm(r_sat)
    
    # Simple cylindrical shadow model
    return sun_sat_dot < 0 and r_sat_mag**2 - sun_sat_dot**2 < R_EARTH**2 