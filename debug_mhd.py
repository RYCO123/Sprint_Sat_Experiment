import numpy as np
from datetime import datetime
from satellite import Satellite
from helpers import calculate_sun_vector, get_geodetic_coords, get_magnetic_field_bgs
from constants import TETHER_LENGTH, TETHER_DIAMETER, TETHER_RESISTIVITY, CIRCUIT_RESISTANCE

def debug_mhd_power():
    """Debug the MHD power calculation to see what's happening."""
    
    # Create MHD satellite
    mhd_sat = Satellite("MHD Sprint", is_sprint_sat=True)
    
    # Get initial state
    state = mhd_sat.initial_conditions
    sim_datetime = datetime.strptime("2024-03-01", "%Y-%m-%d")
    
    print("=== MHD Power Calculation Debug ===")
    print(f"Initial state: {state}")
    print(f"Position (km): {state[:3]/1000}")
    print(f"Velocity (km/s): {state[3:]/1000}")
    
    # Get geodetic coordinates
    lat_deg, lon_deg, alt_m = get_geodetic_coords(state)
    print(f"Latitude: {lat_deg:.2f}°, Longitude: {lon_deg:.2f}°, Altitude: {alt_m/1000:.1f} km")
    
    # Get magnetic field
    b_field_ned_nT = get_magnetic_field_bgs(sim_datetime, lat_deg, lon_deg, alt_m)
    b_ned_T = np.array([b_field_ned_nT['x'], b_field_ned_nT['y'], b_field_ned_nT['z']]) * 1e-9
    print(f"Magnetic field NED (nT): {b_field_ned_nT}")
    print(f"Magnetic field NED (T): {b_ned_T}")
    
    # Transform to ECI coordinates
    lat_rad, lon_rad = np.radians(lat_deg), np.radians(lon_deg)
    clat, slat, clon, slon = np.cos(lat_rad), np.sin(lat_rad), np.cos(lon_rad), np.sin(lon_rad)
    C_ned_to_eci = np.array([[-slat*clon, -slon, -clat*clon], 
                            [-slat*slon, clon, -clat*slon], 
                            [clat, 0.0, -slat]])
    b_eci = C_ned_to_eci @ b_ned_T
    
    # DEBUG: Override with 50 microTesla field pointing downward (nadir)
    r_sat = state[:3]
    r_unit = r_sat / np.linalg.norm(r_sat)
    b_eci = 50e-6 * r_unit  # 50 microTesla pointing toward Earth center
    print(f"Magnetic field ECI (T): {b_eci}")
    
    # Calculate tether vector (nadir-pointing)
    l_tether_vec = -r_sat / np.linalg.norm(r_sat) * TETHER_LENGTH
    print(f"Tether vector (m): {l_tether_vec}")
    
    # Calculate tether resistance
    tether_cross_section = np.pi * (TETHER_DIAMETER / 2)**2
    tether_resistance = TETHER_RESISTIVITY * TETHER_LENGTH / tether_cross_section
    total_resistance = tether_resistance + CIRCUIT_RESISTANCE
    print(f"Tether cross-section: {tether_cross_section:.2e} m²")
    print(f"Tether resistance: {tether_resistance:.2f} Ω")
    print(f"Total resistance: {total_resistance:.2f} Ω")
    
    # Calculate motional EMF
    v_sat = state[3:]
    v_cross_b = np.cross(v_sat, b_eci)
    V_emf = np.dot(v_cross_b, l_tether_vec)
    print(f"Velocity cross B: {v_cross_b}")
    print(f"EMF: {V_emf:.6f} V")
    
    if V_emf > 0:
        I_circuit = V_emf / total_resistance
        power = I_circuit * V_emf
        print(f"Circuit current: {I_circuit:.6f} A")
        print(f"Power: {power:.6f} W")
    else:
        print("EMF is negative or zero - no power generated")
    
    # Test with different positions
    print("\n=== Testing different orbital positions ===")
    test_positions = [
        [6771000, 0, 0],  # 400 km altitude at equator
        [0, 6771000, 0],  # 400 km altitude at 90° longitude
        [0, 0, 6771000],  # 400 km altitude at pole
    ]
    
    for i, pos in enumerate(test_positions):
        test_state = np.array([*pos, 0, 7660, 0])  # Circular orbit velocity
        lat_deg, lon_deg, alt_m = get_geodetic_coords(test_state)
        b_field_ned_nT = get_magnetic_field_bgs(sim_datetime, lat_deg, lon_deg, alt_m)
        b_ned_T = np.array([b_field_ned_nT['x'], b_field_ned_nT['y'], b_field_ned_nT['z']]) * 1e-9
        
        # Transform to ECI
        lat_rad, lon_rad = np.radians(lat_deg), np.radians(lon_deg)
        clat, slat, clon, slon = np.cos(lat_rad), np.sin(lat_rad), np.cos(lon_rad), np.sin(lon_rad)
        C_ned_to_eci = np.array([[-slat*clon, -slon, -clat*clon], 
                                [-slat*slon, clon, -clat*slon], 
                                [clat, 0.0, -slat]])
        b_eci = C_ned_to_eci @ b_ned_T
        
        # DEBUG: Override with 50 microTesla field pointing downward (nadir)
        r_sat = test_state[:3]
        r_unit = r_sat / np.linalg.norm(r_sat)
        b_eci = 50e-6 * r_unit  # 50 microTesla pointing toward Earth center
        
        l_tether_vec = -test_state[:3] / np.linalg.norm(test_state[:3]) * TETHER_LENGTH
        v_cross_b = np.cross(test_state[3:], b_eci)
        V_emf = np.dot(v_cross_b, l_tether_vec)
        
        print(f"Position {i+1}: EMF = {V_emf:.6f} V, Power = {max(0, V_emf**2/total_resistance):.6f} W")

if __name__ == "__main__":
    debug_mhd_power() 