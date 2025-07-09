import numpy as np
from datetime import datetime
from constants import R_EARTH, MU, INITIAL_ALTITUDE, INCLINATION_DEG, DAYS_PER_YEAR, SECONDS_PER_DAY, TWO_PI

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

def get_space_weather_data(date):
    """
    Fetches F10.7 and Ap solar indices from the CelesTrak API.
    
    Args:
        date: datetime object
    
    Returns:
        tuple: (f107_avg, f107_daily, ap_daily)
    """
    try:
        url = "https://celestrak.org/SpaceData/sw19571001.txt"
        import requests
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        lines = response.text.strip().split('\n')
        
        for line in reversed(lines):
            if line.startswith(f"{date.year} {date.month:02d} {date.day:02d}"):
                parts = line.split()
                ap_daily, f107_daily, f107_avg = int(parts[8]), float(parts[12]), float(parts[13])
                result = (f107_avg, f107_daily, ap_daily)
                print(f"Fetched Space Weather for {date.strftime('%Y-%m-%d')}: F10.7_avg={f107_avg}, F10.7_daily={f107_daily}, Ap={ap_daily}")
                return result
    except Exception as e:
        print(f"WARN: Could not fetch space weather data: {e}. Using default values.")
    
    # Default values if API fails
    return (150.0, 150.0, 12)

def calculate_sun_vector(sim_time):
    """
    Calculate sun direction vector (simplified model).
    
    Args:
        sim_time: Simulation time in seconds
    
    Returns:
        np.array: Sun direction vector [x, y, z]
    """
    # Simplified sun position (assume sun is in the x-y plane)
    sun_angle = (sim_time / (DAYS_PER_YEAR * SECONDS_PER_DAY)) * TWO_PI
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