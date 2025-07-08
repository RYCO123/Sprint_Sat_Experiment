import numpy as np
from constants import (MASS_1U, AREA_1U, DRAG_COEFFICIENT, SOLAR_POWER_PEAK, SOLAR_EFFICIENCY)
from helpers import calculate_circular_orbit_conditions, check_eclipse
from ODE import calculate_mhd_generator_power

class Satellite:
    """Satellite class with power generation based on type."""
    
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
        """Calculate power generation based on satellite type."""
        if self.is_sprint_sat:
            return self.calculate_mhd_power(state, sim_datetime, f107a, f107, ap)
        else:
            return self.calculate_solar_power(state, sun_vector)
    
    def calculate_solar_power(self, state, sun_vector):
        """Calculate solar power generation with eclipse detection."""
        if check_eclipse(state, sun_vector):
            return 0.0  # No power in eclipse
        
        # Calculate power based on sun angle
        r_sat = state[:3]
        panel_normal = -r_sat / np.linalg.norm(r_sat)  # Nadir-pointing panels
        cos_theta = np.dot(panel_normal, sun_vector)
        
        return SOLAR_POWER_PEAK * SOLAR_EFFICIENCY * max(0, cos_theta)
    
    def calculate_mhd_power(self, state, sim_datetime, f107a, f107, ap):
        """Calculate MHD power generation using MHD generator."""
        return calculate_mhd_generator_power(state, sim_datetime, f107a, f107, ap)
    
    def log_step(self, time, state, power):
        """Log simulation step data."""
        self.time_history.append(time)
        self.state_history.append(state)
        
        # Ensure power is a scalar
        power_scalar = float(power) if hasattr(power, '__iter__') else float(power)
        self.power_history.append(power_scalar)
        
        # Update total energy
        if len(self.time_history) > 1:
            dt = time - self.time_history[-2]
            self.total_energy_J += power_scalar * dt
        
        self.energy_history.append(self.total_energy_J)
    
    def get_current_altitude(self):
        """Get current altitude from latest state."""
        if not self.state_history:
            return 0.0
        
        current_state = self.state_history[-1]
        return np.linalg.norm(current_state[:3]) - 6371000  # R_EARTH
    
    def get_average_power(self):
        """Calculate average power over the mission."""
        return np.mean(self.power_history) if self.power_history else 0.0
    
    def get_total_energy(self):
        """Get total energy generated."""
        return self.total_energy_J 