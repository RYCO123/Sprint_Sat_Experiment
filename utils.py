"""
Utility functions for Sprint Satellite simulation.

This module contains helper functions for space weather data,
orbital calculations, and other utility operations.
"""

import numpy as np
from datetime import datetime
import requests
from config import R_EARTH, INITIAL_ALTITUDE, INCLINATION_DEG, MU, TWO_PI, DAYS_PER_YEAR, SECONDS_PER_DAY


def get_space_weather_data(date):
    """
    Fetch F10.7 and Ap solar indices from the CelesTrak API.
    
    Args:
        date: datetime object for the requested date
        
    Returns:
        tuple: (f107_avg, f107_daily, ap_daily)
    """
    try:
        url = "https://celestrak.org/SpaceData/sw19571001.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        lines = response.text.strip().split('\n')
        
        for line in reversed(lines):
            if line.startswith(f"{date.year} {date.month:02d} {date.day:02d}"):
                parts = line.split()
                ap_daily = int(parts[8])
                f107_daily = float(parts[12])
                f107_avg = float(parts[13])
                return (f107_avg, f107_daily, ap_daily)
        
        # If no data found, return default values
        return (150.0, 150.0, 12)
        
    except Exception as e:
        print(f"Could not fetch space weather data: {e}. Using default values.")
        return (150.0, 150.0, 12)


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