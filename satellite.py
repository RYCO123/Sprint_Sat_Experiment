import numpy as np
from constants import (MASS_1U, AREA_1U, DRAG_COEFFICIENT, TETHER_LENGTH, TETHER_DIAMETER,
                      TETHER_RESISTIVITY, CIRCUIT_RESISTANCE, SOLAR_POWER_PEAK, SOLAR_EFFICIENCY)
from helpers import (calculate_circular_orbit_conditions, get_geodetic_coords, 
                    get_magnetic_field_bgs, check_eclipse)

class Satellite:
    """
    Satellite class with simplified power generation based on type.
    
    Args:
        name (str): Satellite name
        is_sprint_sat (bool): True for MHD sprint satellite, False for standard solar
    """
    
    def __init__(self, name, is_sprint_sat=False):
        self.name = name
        self.is_sprint_sat = is_sprint_sat
        self.mass = MASS_1U
        self.area = AREA_1U
        self.c_d = DRAG_COEFFICIENT
        
        # Initialize orbital state
        self.initial_conditions = calculate_circular_orbit_conditions()
        
        # Data storage
        self.state_history = []
        self.time_history = []
        self.power_history = []
        self.energy_history = []
        self.total_energy_J = 0.0
    
    def calculate_power(self, state, sim_datetime, sun_vector, f107a, f107, ap):
        """
        Calculate power generation based on satellite type.
        
        Args:
            state: State vector [x, y, z, vx, vy, vz]
            sim_datetime: Simulation datetime
            sun_vector: Sun direction vector
            f107a, f107, ap: Space weather parameters
        
        Returns:
            float: Power generated in Watts
        """
        if self.is_sprint_sat:
            return self.calculate_mhd_power(state, sim_datetime, f107a, f107, ap)
        else:
            return self.calculate_solar_power(state, sun_vector)
    
    def calculate_solar_power(self, state, sun_vector):
        """
        Calculate solar power generation with eclipse detection.
        
        Args:
            state: State vector [x, y, z, vx, vy, vz]
            sun_vector: Sun direction vector
        
        Returns:
            float: Solar power in Watts
        """
        # Check if satellite is in eclipse
        if check_eclipse(state, sun_vector):
            return 0.0  # No power in eclipse
        
        # Calculate power based on sun angle
        r_sat = state[:3]
        panel_normal = -r_sat / np.linalg.norm(r_sat)  # Nadir-pointing panels
        cos_theta = np.dot(panel_normal, sun_vector)
        
        # Power calculation with efficiency
        power = SOLAR_POWER_PEAK * SOLAR_EFFICIENCY * max(0, cos_theta)
        
        return power
    
    def calculate_mhd_power(self, state, sim_datetime, f107a, f107, ap):
        """
        Calculate MHD power generation using deployable electrodynamic tether.
        
        Args:
            state: State vector [x, y, z, vx, vy, vz]
            sim_datetime: Simulation datetime
            f107a, f107, ap: Space weather parameters
        
        Returns:
            float: MHD power in Watts
        """
        lat_deg, lon_deg, alt_m = get_geodetic_coords(state)
        

        b_field_ned_nT = get_magnetic_field_bgs(sim_datetime, lat_deg, lon_deg, alt_m)
        b_ned_T = np.array([b_field_ned_nT['x'], b_field_ned_nT['y'], b_field_ned_nT['z']]) * 1e-9
        
        # Transform to ECI coordinates
        lat_rad, lon_rad = np.radians(lat_deg), np.radians(lon_deg)
        clat, slat, clon, slon = np.cos(lat_rad), np.sin(lat_rad), np.cos(lon_rad), np.sin(lon_rad)
        C_ned_to_eci = np.array([[-slat*clon, -slon, -clat*clon], 
                                [-slat*slon, clon, -clat*slon], 
                                [clat, 0.0, -slat]])
        b_eci = C_ned_to_eci @ b_ned_T
        

        r_sat = state[:3]
        r_unit = r_sat / np.linalg.norm(r_sat)
        b_eci = 50e-6 * r_unit
        
        # Calculate tether vector (perpendicular to velocity and magnetic field for power generation)
        v_sat = state[3:]
        
        # Make tether perpendicular to both velocity and magnetic field
        # This ensures maximum power generation
        tether_direction = np.cross(v_sat, b_eci)
        if np.linalg.norm(tether_direction) > 0:
            tether_direction = tether_direction / np.linalg.norm(tether_direction)
        else:
            # Fallback to nadir direction if cross product is zero
            tether_direction = -r_sat / np.linalg.norm(r_sat)
        
        l_tether_vec = tether_direction * TETHER_LENGTH
        
        # Calculate tether resistance based on physical properties
        tether_cross_section = np.pi * (TETHER_DIAMETER / 2)**2  # m²
        tether_resistance = TETHER_RESISTIVITY * TETHER_LENGTH / tether_cross_section  # Ohms
        total_resistance = tether_resistance + CIRCUIT_RESISTANCE  # Total circuit resistance
        
        # Calculate motional EMF
        v_cross_b = np.cross(v_sat, b_eci)
        V_emf = np.dot(v_cross_b, l_tether_vec)
        
        if V_emf > 0:
            I_circuit = V_emf / total_resistance
            power = I_circuit * V_emf
            return power
        else:
            return 0.0
    
    def log_step(self, time, state, power):
        """
        Log simulation step data.
        
        Args:
            time: Current simulation time
            state: Current state vector
            power: Current power generation
        """
        self.time_history.append(time)
        self.state_history.append(state)
        self.power_history.append(power)
        
        # Update total energy
        if len(self.time_history) > 1:
            dt = time - self.time_history[-2]
            self.total_energy_J += power * dt
        
        self.energy_history.append(self.total_energy_J)
    
    def get_current_altitude(self):
        """
        Get current altitude from latest state.
        
        Returns:
            float: Current altitude in meters
        """
        if not self.state_history:
            return 0.0
        
        current_state = self.state_history[-1]
        return np.linalg.norm(current_state[:3]) - 6371000  # R_EARTH
    
    def get_average_power(self):
        """
        Calculate average power over the mission.
        
        Returns:
            float: Average power in Watts
        """
        if not self.power_history:
            return 0.0
        
        return np.mean(self.power_history)
    
    def get_total_energy(self):
        """
        Get total energy generated.
        
        Returns:
            float: Total energy in Joules
        """
        return self.total_energy_J 