import numpy as np
from constants import MU, R_EARTH, TETHER_LENGTH, TETHER_DIAMETER, TETHER_RESISTIVITY, CIRCUIT_RESISTANCE
from helpers import get_geodetic_coords, get_magnetic_field_bgs
from astropy.coordinates import EarthLocation, ITRS, GCRS, TETE
from astropy.time import Time
import astropy.units as u

def orbital_ode_standard(t, state, satellite, sim_datetime, sun_vector, f107a, f107, ap):
    """
    Orbital dynamics ODE for standard solar satellite (no drag terms).
    
    Args:
        t: Time (s)
        state: State vector [x, y, z, vx, vy, vz]
        satellite: Satellite object
        sim_datetime: Simulation datetime
        sun_vector: Sun direction vector
        f107a, f107, ap: Space weather parameters (unused for standard sat)
    
    Returns:
        Derivatives [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]
    """
    # Extract position and velocity
    r = state[:3]
    v = state[3:]
    
    # Gravitational acceleration only (no drag for standard satellite)
    r_mag = np.linalg.norm(r)
    a_grav = -MU * r / (r_mag**3)
    
    # Return derivatives [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]
    return np.concatenate([v, a_grav])

def orbital_ode_sprint(t, state, satellite, sim_datetime, sun_vector, f107a, f107, ap):
    """
    Orbital dynamics ODE for MHD sprint satellite (includes MHD drag calculation).
    
    Args:
        t: Time (s)
        state: State vector [x, y, z, vx, vy, vz]
        satellite: Satellite object
        sim_datetime: Simulation datetime
        sun_vector: Sun direction vector
        f107a, f107, ap: Space weather parameters for magnetic field calculation
    
    Returns:
        Derivatives [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]
    """
    # Extract position and velocity
    r = state[:3]
    v = state[3:]
    
    # Gravitational acceleration
    r_mag = np.linalg.norm(r)
    a_grav = -MU * r / (r_mag**3)
    
    # Initialize total acceleration
    a_total = a_grav
    
    # Calculate MHD drag force directly in the ODE
    f_mhd = calculate_mhd_drag_in_ode(state, sim_datetime, f107a, f107, ap)
    a_mhd = f_mhd / satellite.mass
    a_total += a_mhd
    
    # Return derivatives [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]
    return np.concatenate([v, a_total])

def calculate_mhd_drag_in_ode(state, sim_datetime, f107a, f107, ap):
    """
    Calculate MHD drag force vector directly in the ODE using deployable tether.
    
    Args:
        state: State vector [x, y, z, vx, vy, vz]
        sim_datetime: Simulation datetime
        f107a, f107, ap: Space weather parameters
    
    Returns:
        np.array: MHD drag force vector [Fx, Fy, Fz]
    """
    lat_deg, lon_deg, alt_m = get_geodetic_coords(state)
    
    # Get magnetic field
    b_field_ned_nT = get_magnetic_field_bgs(sim_datetime, lat_deg, lon_deg, alt_m)
    b_ned_T = np.array([b_field_ned_nT['x'], b_field_ned_nT['y'], b_field_ned_nT['z']]) * 1e-9
    
    # Transform to ECI coordinates
    lat_rad, lon_rad = np.radians(lat_deg), np.radians(lon_deg)




    clat, slat, clon, slon = np.cos(lat_rad), np.sin(lat_rad), np.cos(lon_rad), np.sin(lon_rad)
    C_ned_to_eci = np.array([[-slat*clon, -slon, -clat*clon], 
                            [-slat*slon, clon, -clat*slon], 
                            [clat, 0.0, -slat]])

    b_ecef = C_ned_to_eci @ b_ned_T * u.T # Add units for astropy

    location = EarthLocation.from_geodetic(lon_deg*u.deg, lat_deg*u.deg, alt_m*u.m)
    time = Time(sim_datetime)

    frame_itrs = ITRS(x=location.x, y=location.y, z=location.z, 
                    v_x=b_ecef[0], v_y=b_ecef[1], v_z=b_ecef[2],
                    obstime=time)

    frame_gcrs = frame_itrs.transform_to(GCRS(obstime=time))

    b_eci = frame_gcrs.velocity.d_xyz.to(u.T).value                   

    
    # Calculate tether vector (nadir-pointing)
    r_sat = state[:3]
    l_tether_vec = -r_sat / np.linalg.norm(r_sat) * TETHER_LENGTH
    
    # Calculate tether resistance based on physical properties
    tether_cross_section = np.pi * (TETHER_DIAMETER / 2)**2  # m²
    tether_resistance = TETHER_RESISTIVITY * TETHER_LENGTH / tether_cross_section  # Ohms
    total_resistance = tether_resistance + CIRCUIT_RESISTANCE  # Total circuit resistance
    
    # Calculate motional EMF and current
    v_sat = state[3:]
    v_cross_b = np.cross(v_sat, b_eci)
    V_emf = np.dot(v_cross_b, l_tether_vec)
    
    if V_emf > 0:
        I_circuit = V_emf / total_resistance
        # MHD drag force: F = I * (L × B)
        f_mhd = I_circuit * np.cross(l_tether_vec, b_eci)
        return f_mhd
    else:
        return np.zeros(3)

def get_ode_function(is_sprint_sat):
    """
    Get the appropriate ODE function based on satellite type.
    
    Args:
        is_sprint_sat: Boolean indicating if this is a sprint satellite
    
    Returns:
        ODE function (orbital_ode_standard or orbital_ode_sprint)
    """
    if is_sprint_sat:
        return orbital_ode_sprint
    else:
        return orbital_ode_standard