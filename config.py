"""
Configuration file for Sprint Satellite simulation.

This module contains all physical constants, simulation parameters,
and configuration values used throughout the simulation.
"""

import numpy as np

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

# Gravitational constants
G = 6.67430e-11  # Gravitational constant (m³/kg/s²)
M_EARTH = 5.972e24  # Earth's mass (kg)
R_EARTH = 6371000  # Earth's radius (m)
MU = G * M_EARTH  # Earth's gravitational parameter (m³/s²)

# Plasma physics constants
ELECTRON_CHARGE = 1.602e-19  # C
ELECTRON_MASS = 9.109e-31  # kg
BOLTZMANN_CONSTANT = 1.381e-23  # J/K
PLASMA_TEMPERATURE = 1000  # K (typical ionospheric temperature)

# Time constants
SECONDS_PER_DAY = 86400.0  # seconds in a day
DAYS_PER_YEAR = 365.25  # days in a year (including leap years)
TWO_PI = 2 * np.pi  # 2π for angular calculations

# =============================================================================
# SATELLITE PARAMETERS
# =============================================================================

# 1U CubeSat parameters
MASS_1U = 1.33  # kg

# Solar satellite parameters
SOLAR_POWER_PEAK = 10.0  # W (peak power)
SOLAR_EFFICIENCY = 0.25  # 25% average efficiency

# MHD Generator parameters
ELECTRODE_DISTANCE = 0.1  # m (distance between perpendicular electrodes)
GENERATOR_LENGTH = 0.1135  # m (longitudinal distance of the generator)
MAGNETIC_DISTANCE = 0.05  # m (distance between magnets creating B field)
MHD_MAGNET_STRENGTH = 0.0028  # Tesla (magnetic field strength from magnets)
MHD_CONDUCTOR_RESISTIVITY = 1.68e-8  # Ω⋅m (copper resistivity)
MHD_CONDUCTOR_DIAMETER = 0.001  # m (1 mm diameter conductors)
MHD_CIRCUIT_RESISTANCE = 5.0  # Ohms (additional circuit resistance)

# =============================================================================
# ORBITAL PARAMETERS
# =============================================================================

INITIAL_ALTITUDE = 400.0e3  # m (400 km)
INCLINATION_DEG = 20.0  # degrees
TERMINATION_ALTITUDE = 60.1e3  # m (60.1 km - deorbit boundary)

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

SIMULATION_START_DATE = "2024-03-01"

# Space weather parameters (default values)
F107A = 150.0  # F10.7 81-day average (solar flux units)
F107 = 150.0   # F10.7 daily (solar flux units)
AP = 4.0       # Ap index (geomagnetic activity) 