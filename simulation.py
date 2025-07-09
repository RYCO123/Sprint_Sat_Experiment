"""
Core simulation logic for Sprint Satellite simulation.

This module contains the main Simulation class that orchestrates
the orbital dynamics, power calculations, and data processing.
"""

import numpy as np
from scipy.integrate import solve_ivp
from datetime import datetime
from tqdm.auto import tqdm
import os

from config import SIMULATION_START_DATE
from satellite import Satellite
from physics import get_ode_function
from utils import get_space_weather_data, calculate_sun_vector
from plotting import (plot_orbital_trajectories_2d, plot_orbital_trajectories_3d,
                     create_comprehensive_comparison)
from analysis import generate_summary, predict_deorbit_time


class Simulation:
    """
    Main simulation class that runs both satellites and generates comparisons.
    """
    
    def __init__(self):
        """Initialize the simulation with satellites and data storage."""
        # Create satellites
        self.mhd_sat = Satellite("MHD Sprint Satellite", is_sprint_sat=True)
        self.solar_sat = Satellite("Standard Solar Satellite", is_sprint_sat=False)
        
        # Simulation data storage
        self.time_data = []
        self.mhd_states = []
        self.solar_states = []
        self.mhd_power = []
        self.solar_power = []
        self.mhd_energy = []
        self.solar_energy = []
    
    def run_simulation(self, t_span, num_steps=150):
        """
        Run the orbital simulation using ODE solver with discrete time steps.
        
        Args:
            t_span: Tuple of (start_time, end_time) in seconds
            num_steps: Number of discrete time steps (default: 150)
        """
        print("Starting Sprint Satellite vs Standard Solar Satellite simulation...")
        print(f"Simulation duration: {t_span[1]/3600:.1f} hours")
        print(f"Number of time steps: {num_steps}")
        
        # Get initial space weather data
        sim_datetime = datetime.strptime(SIMULATION_START_DATE, "%Y-%m-%d")
        f107a, f107, ap = get_space_weather_data(sim_datetime)
        
        # Use direct PyIRI calls for plasma density and magnetic field
        print("Using direct PyIRI calls for ionospheric data")
        
        # Get ODE functions
        mhd_ode = get_ode_function(is_sprint_sat=True)
        solar_ode = get_ode_function(is_sprint_sat=False)
        
        # Initial conditions
        y0_mhd = self.mhd_sat.initial_conditions
        y0_solar = self.solar_sat.initial_conditions
        
        # Create discrete time array
        t_eval = np.linspace(t_span[0], t_span[1], num_steps)
        
        print("Solving ODE for MHD Sprint Satellite...")
        # Solve ODE for MHD satellite
        with tqdm(total=t_span[1], desc="MHD Sprint Sat", unit="s") as pbar:
            # This wrapper function will update the progress bar
            def mhd_ode_with_progress(t, y, *args):
                # Update the progress bar to the current time `t`
                pbar.update(t - pbar.n)
                # Call the original ODE function
                return mhd_ode(t, y, *args)

            # Solve ODE for MHD satellite, using our new wrapper function
            solution_mhd = solve_ivp(
                fun=lambda t, y: mhd_ode_with_progress(t, y, self.mhd_sat, sim_datetime, 
                                        calculate_sun_vector(t), f107a, f107, ap),
                t_span=t_span,
                y0=y0_mhd,
                method='RK45',
                t_eval=t_eval,
                rtol=1e-8,
                atol=1e-10
            )
        
        print("Solving ODE for Standard Solar Satellite...")
        with tqdm(total=t_span[1], desc="Standard Solar Sat", unit="s") as pbar:
            # This wrapper function will update the progress bar
            def solar_ode_with_progress(t, y, *args):
                # Update the progress bar to the current time `t`
                pbar.update(t - pbar.n)
                # Call the original ODE function
                return solar_ode(t, y, *args)

            # Solve ODE for solar satellite, using our new wrapper function
            solution_solar = solve_ivp(
                fun=lambda t, y: solar_ode_with_progress(t, y, self.solar_sat, sim_datetime, 
                                        calculate_sun_vector(t), f107a, f107, ap),
                t_span=t_span,
                y0=y0_solar,
                method='RK45',
                t_eval=t_eval,
                rtol=1e-8,
                atol=1e-10
            )
        
        # Process results
        self.process_results(solution_mhd, solution_solar, sim_datetime, f107a, f107, ap)
        
        print("Simulation completed successfully!")
    
    def process_results(self, solution_mhd, solution_solar, sim_datetime, f107a, f107, ap):
        """
        Process simulation results and calculate power/energy.
        
        Args:
            solution_mhd: ODE solution for MHD satellite
            solution_solar: ODE solution for solar satellite
            sim_datetime: Simulation datetime
            f107a, f107, ap: Space weather parameters
        """
        times = solution_mhd.t
        
        print("Processing results and calculating power generation...")
        
        for i, t in enumerate(tqdm(times, desc="Processing simulation data")):
            # Get states
            state_mhd = solution_mhd.y[:, i]
            state_solar = solution_solar.y[:, i]
            
            # Calculate sun vector
            sun_vector = calculate_sun_vector(t)
            
            # Calculate powers
            power_mhd = self.mhd_sat.calculate_power(state_mhd, sim_datetime, sun_vector, f107a, f107, ap)
            power_solar = self.solar_sat.calculate_power(state_solar, sim_datetime, sun_vector, f107a, f107, ap)
            
            # Ensure powers are scalars
            power_mhd_scalar = float(power_mhd) if hasattr(power_mhd, '__iter__') else float(power_mhd)
            power_solar_scalar = float(power_solar) if hasattr(power_solar, '__iter__') else float(power_solar)
            
            # Log data
            self.mhd_sat.log_step(t, state_mhd, power_mhd_scalar)
            self.solar_sat.log_step(t, state_solar, power_solar_scalar)
            
            # Store data for plotting
            self.time_data.append(t)
            self.mhd_states.append(state_mhd)
            self.solar_states.append(state_solar)
            self.mhd_power.append(power_mhd_scalar)
            self.solar_power.append(power_solar_scalar)
            self.mhd_energy.append(float(self.mhd_sat.total_energy_J))
            self.solar_energy.append(float(self.solar_sat.total_energy_J))
    
    def generate_plots(self):
        """Generate comparison plots."""
        print("Generating comparison plots...")

        if not os.path.exists('plots'):
            os.makedirs('plots')

        # 2D orbital trajectories
        plot_orbital_trajectories_2d(self.mhd_states, self.solar_states, 'plots/orbital_trajectories_2d.png')
        
        # 3D orbital trajectories
        plot_orbital_trajectories_3d(self.mhd_states, self.solar_states, 'plots/orbital_trajectories_3d.png')
        
        # Comprehensive comparison
        create_comprehensive_comparison(
            self.time_data, self.mhd_states, self.solar_states,
            self.mhd_power, self.solar_power, self.mhd_energy, self.solar_energy,
            'plots/comprehensive_comparison.png'
        )
        
        print("All plots generated successfully!")
    
    def generate_summary(self):
        """Generate simulation summary and analysis."""
        generate_summary(
            self.time_data, self.mhd_states, self.solar_states,
            self.mhd_power, self.solar_power, self.mhd_energy, self.solar_energy
        )
        
        # Perform polynomial regression to predict deorbit time
        predict_deorbit_time(self.time_data, self.mhd_states) 