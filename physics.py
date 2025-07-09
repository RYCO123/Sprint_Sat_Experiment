"""
Physics calculations for Sprint Satellite simulation.

This module contains all physics calculations including orbital dynamics,
MHD generator physics, plasma physics, and ionospheric modeling.
"""

import numpy as np
import PyIRI
import PyIRI.edp_update as ml
from config import (MU, R_EARTH, ELECTRODE_DISTANCE, GENERATOR_LENGTH, MAGNETIC_DISTANCE,
                   MHD_MAGNET_STRENGTH, MHD_CONDUCTOR_RESISTIVITY, 
                   MHD_CONDUCTOR_DIAMETER, MHD_CIRCUIT_RESISTANCE,
                   ELECTRON_CHARGE, ELECTRON_MASS, BOLTZMANN_CONSTANT, PLASMA_TEMPERATURE)


def orbital_ode_standard(t, state, satellite, sim_datetime, sun_vector, f107a, f107, ap):
    """
    Standard solar satellite ODE - gravitational forces only.
    
    Args:
        t: Time (unused but required by ODE solver)
        state: State vector [x, y, z, vx, vy, vz]
        satellite: Satellite object
        sim_datetime: Simulation datetime
        sun_vector: Sun direction vector
        f107a, f107, ap: Space weather parameters (unused)
        
    Returns:
        np.ndarray: State derivative [vx, vy, vz, ax, ay, az]
    """
    r = state[:3]
    v = state[3:]
    
    # Gravitational acceleration only
    r_mag = np.linalg.norm(r)
    a_grav = -MU * r / (r_mag**3)
    
    return np.concatenate([v, a_grav])


def orbital_ode_sprint(t, state, satellite, sim_datetime, sun_vector, f107a, f107, ap):
    """
    MHD sprint satellite ODE - includes MHD generator drag forces.
    
    Args:
        t: Time (unused but required by ODE solver)
        state: State vector [x, y, z, vx, vy, vz]
        satellite: Satellite object
        sim_datetime: Simulation datetime
        sun_vector: Sun direction vector
        f107a, f107, ap: Space weather parameters
        
    Returns:
        np.ndarray: State derivative [vx, vy, vz, ax, ay, az]
    """
    r = state[:3]
    v = state[3:]
    
    # Gravitational acceleration
    r_mag = np.linalg.norm(r)
    a_grav = -MU * r / (r_mag**3)
    
    # MHD generator drag force
    f_mhd = calculate_mhd_generator_drag_force(state, sim_datetime, f107a, f107, ap)
    a_mhd = f_mhd / satellite.mass
    
    return np.concatenate([v, a_grav + a_mhd])


def get_ionospheric_data(state, sim_datetime, f107a, f107, ap):
    """
    Get ionospheric plasma density and magnetic field using PyIRI.
    
    Args:
        state: State vector [x, y, z, vx, vy, vz]
        sim_datetime: Simulation datetime
        f107a, f107, ap: Space weather parameters
        
    Returns:
        tuple: (electron_density_m3, magnetic_field_T)
    """
    r = state[:3]
    altitude_km = (np.linalg.norm(r) - R_EARTH) / 1000.0

    # Convert position to latitude/longitude
    lat_rad = np.arcsin(r[2] / np.linalg.norm(r))
    lon_rad = np.arctan2(r[1], r[0])
    lat_deg = np.degrees(lat_rad)
    lon_deg = np.degrees(lon_rad)

    # IRI valid range: 60-2000 km
    altitude_km = np.clip(altitude_km, 60.0, 2000.0)

    alon = np.array([lon_deg])  # Single location longitude
    alat = np.array([lat_deg])  # Single location latitude
    aalt = np.array([altitude_km])  # Single altitude point
    ahr = np.array([sim_datetime.hour + sim_datetime.minute/60 + sim_datetime.second/3600])
    
    year = sim_datetime.year
    month = sim_datetime.month
    day = sim_datetime.day
    ccir_or_ursi = 0  # 0=CCIR, 1=URSI

    f2, f1, e_peak, es_peak, sun, mag, edp = ml.IRI_density_1day(
        year, month, day, ahr, alon, alat, aalt, f107, PyIRI.coeff_dir, ccir_or_ursi
    )
    
    # Extract electron density from edp array
    n_e = edp[0, 0, 0]
    n_e = max(n_e, 1e9)  # Minimum density
    
    inc = mag['inc'][0]  # Magnetic field inclination in degrees
    modip = mag['modip'][0] # Modified dip angle in degrees
    mag_dip_lat = mag['mag_dip_lat'][0]  # Magnetic dip latitude in degrees
    
    inc_rad = np.radians(inc)
    modip_rad = np.radians(modip)
    
    r_mag = np.linalg.norm(r)
    B_magnitude = 3.12e-5 * (R_EARTH / r_mag)**3  # Approximate Earth's dipole field
    
    B_horizontal = B_magnitude * np.cos(inc_rad)
    B_vertical = B_magnitude * np.sin(inc_rad)
    
    B_T = np.array([0, 0, -B_magnitude])  # Simplified southward field
    
    return n_e, B_T


def calculate_plasma_conductivity(n_e, T_e):
    """
    Calculate plasma conductivity using Spitzer formula.
    
    Args:
        n_e: Electron density (m^-3)
        T_e: Electron temperature (K)
        
    Returns:
        float: Plasma conductivity (S/m)
    """
    ln_Lambda = 10.0
    epsilon_0 = 8.854e-12  # F/m
    nu_ei = (4 * np.pi * n_e * ELECTRON_CHARGE**4 * ln_Lambda) / \
            (4 * np.pi * epsilon_0**2 * np.sqrt(ELECTRON_MASS) * (BOLTZMANN_CONSTANT * T_e)**1.5)
    
    tau = 1.0 / nu_ei if nu_ei > 0 else 1e-6
    sigma = n_e * ELECTRON_CHARGE**2 * tau / ELECTRON_MASS
    
    return sigma


def calculate_mhd_generator_drag_force(state, sim_datetime, f107a, f107, ap):
    """
    Calculate MHD generator drag force.
    
    Args:
        state: State vector [x, y, z, vx, vy, vz]
        sim_datetime: Simulation datetime
        f107a, f107, ap: Space weather parameters
        
    Returns:
        np.ndarray: Drag force vector (N)
    """
    # Get ionospheric data
    n_e, B_earth = get_ionospheric_data(state, sim_datetime, f107a, f107, ap)
    
    # MHD generator configuration
    B_magnets = np.array([0, 0, MHD_MAGNET_STRENGTH])
    B_total = B_earth + B_magnets
    
    # Plasma properties
    sigma_plasma = calculate_plasma_conductivity(n_e, PLASMA_TEMPERATURE)
    
    # Electrode orientation
    velocity_vector = state[3:]
    v_mag = np.linalg.norm(velocity_vector)
    
    if v_mag > 0:
        v_unit = velocity_vector / v_mag
        B_unit = B_total / np.linalg.norm(B_total)
        electrode_orientation = np.cross(v_unit, B_unit)
        electrode_mag = np.linalg.norm(electrode_orientation)
        
        if electrode_mag > 0:
            electrode_orientation = electrode_orientation / electrode_mag
        else:
            electrode_orientation = np.array([1, 0, 0])
    else:
        electrode_orientation = np.array([1, 0, 0])
    
    # Motional EMF calculation
    L_vector = ELECTRODE_DISTANCE * electrode_orientation
    motional_E_field = np.cross(velocity_vector, B_total)
    V_emf = np.dot(motional_E_field, L_vector)
    
    # System resistance calculation
    conductor_cross_section = np.pi * (MHD_CONDUCTOR_DIAMETER / 2)**2
    total_conductor_length = 2 * GENERATOR_LENGTH + 2 * MAGNETIC_DISTANCE
    R_conductor = MHD_CONDUCTOR_RESISTIVITY * total_conductor_length / conductor_cross_section
    
    plasma_area = GENERATOR_LENGTH * MAGNETIC_DISTANCE
    R_plasma = ELECTRODE_DISTANCE / (sigma_plasma * plasma_area)
    R_total = R_conductor + R_plasma + MHD_CIRCUIT_RESISTANCE
    
    # Current and drag force
    if V_emf > 0:
        I = V_emf / R_total
        F_drag = I * np.cross(L_vector, B_total)
    else:
        F_drag = np.zeros(3)
    
    return F_drag


def calculate_mhd_generator_power(state, sim_datetime, f107a, f107, ap):
    """
    Calculate power generated by MHD generator.
    
    Args:
        state: State vector [x, y, z, vx, vy, vz]
        sim_datetime: Simulation datetime
        f107a, f107, ap: Space weather parameters
        
    Returns:
        float: Power generated (W)
    """
    # Get ionospheric data
    n_e, B_earth = get_ionospheric_data(state, sim_datetime, f107a, f107, ap)
    
    # MHD generator configuration
    B_magnets = np.array([0, 0, MHD_MAGNET_STRENGTH])
    B_total = B_earth + B_magnets
    sigma_plasma = calculate_plasma_conductivity(n_e, PLASMA_TEMPERATURE)
    
    # Electrode orientation
    velocity_vector = state[3:]
    v_mag = np.linalg.norm(velocity_vector)
    
    if v_mag > 0:
        v_unit = velocity_vector / v_mag
        B_unit = B_total / np.linalg.norm(B_total)
        electrode_orientation = np.cross(v_unit, B_unit)
        electrode_mag = np.linalg.norm(electrode_orientation)
        
        if electrode_mag > 0:
            electrode_orientation = electrode_orientation / electrode_mag
        else:
            electrode_orientation = np.array([1, 0, 0])
    else:
        electrode_orientation = np.array([1, 0, 0])
    
    # Power calculation
    L_vector = ELECTRODE_DISTANCE * electrode_orientation
    motional_E_field = np.cross(velocity_vector, B_total)
    V_emf = np.dot(motional_E_field, L_vector)
    
    # System resistance
    conductor_cross_section = np.pi * (MHD_CONDUCTOR_DIAMETER / 2)**2
    total_conductor_length = 2 * GENERATOR_LENGTH + 2 * MAGNETIC_DISTANCE
    R_conductor = MHD_CONDUCTOR_RESISTIVITY * total_conductor_length / conductor_cross_section
    plasma_area = GENERATOR_LENGTH * MAGNETIC_DISTANCE
    R_plasma = ELECTRODE_DISTANCE / (sigma_plasma * plasma_area)
    R_total = R_conductor + R_plasma + MHD_CIRCUIT_RESISTANCE
    
    # Power generation
    if V_emf > 0:
        I = V_emf / R_total
        P_generated = V_emf * I
    else:
        P_generated = 0.0
    
    return P_generated


def get_ode_function(is_sprint_sat):
    """
    Get appropriate ODE function based on satellite type.
    
    Args:
        is_sprint_sat: True if MHD Sprint Satellite, False if standard solar
        
    Returns:
        function: ODE function to use
    """
    return orbital_ode_sprint if is_sprint_sat else orbital_ode_standard 