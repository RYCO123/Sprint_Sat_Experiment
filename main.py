import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from datetime import datetime, timedelta
from tqdm.auto import tqdm
from datetime import datetime
import os
from constants import TIME_STEP, SIMULATION_START_DATE, TERMINATION_ALTITUDE
from satellite import Satellite
from ODE import get_ode_function
from helpers import get_space_weather_data, calculate_sun_vector
from plotting import (plot_power_comparison, plot_total_energy_comparison, 
                     plot_orbital_trajectories_2d, plot_orbital_trajectories_3d,
                     plot_altitude_comparison, create_comprehensive_comparison)

class Simulation:
    """
    Main simulation class that runs both satellites and generates comparisons.
    """
    
    def __init__(self):
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
        """Generate simulation summary and save to CSV."""
        print("\n--- Simulation Summary ---")
        
        if not self.time_data:
            print("No simulation data available.")
            return
        
        # Calculate final metrics
        final_time = self.time_data[-1]
        mission_duration_hours = final_time / 3600
        
        # Final altitudes
        final_mhd_state = self.mhd_states[-1]
        final_solar_state = self.solar_states[-1]
        final_mhd_alt = np.linalg.norm(final_mhd_state[:3]) - 6371000
        final_solar_alt = np.linalg.norm(final_solar_state[:3]) - 6371000
        
        # Average powers
        avg_mhd_power = np.mean(self.mhd_power)
        avg_solar_power = np.mean(self.solar_power)
        
        # Total energies
        total_mhd_energy = self.mhd_energy[-1] / 3600  # Convert to Wh
        total_solar_energy = self.solar_energy[-1] / 3600  # Convert to Wh
        
        summary_data = {
            'Metric': [
                'Mission Duration (hours)',
                'Initial Altitude (km)',
                'Final Altitude (km)',
                'Average Power Generated (W)',
                'Total Energy Generated (Wh)',
                'Power Ratio (MHD/Solar)',
                'Energy Ratio (MHD/Solar)'
            ],
            'MHD Sprint Satellite': [
                f"{mission_duration_hours:.5f}",
                f"400.00000",
                f"{final_mhd_alt/1000:.5f}",
                f"{avg_mhd_power:.5f}",
                f"{total_mhd_energy:.5f}",
                f"{avg_mhd_power/avg_solar_power:.5f}x",
                f"{total_mhd_energy/total_solar_energy:.5f}x"
            ],
            'Standard Solar Satellite': [
                f"{mission_duration_hours:.5f}",
                f"400.00000",
                f"{final_solar_alt/1000:.5f}",
                f"{avg_solar_power:.5f}",
                f"{total_solar_energy:.5f}",
                "1.00000x",
                "1.00000x"
            ]
        }
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv('simulation_summary.csv', index=False)
        print(summary_df.to_string(index=False))
        print("\nSummary saved to simulation_summary.csv")
        
        # Print key insights
        print(f"\n--- Key Insights ---")
        print(f"MHD Sprint Satellite generates {avg_mhd_power/avg_solar_power:.1f}x more power")
        print(f"MHD Sprint Satellite generates {total_mhd_energy/total_solar_energy:.1f}x more total energy")
        print(f"Simulation completed in {mission_duration_hours:.1f} hours with {len(self.time_data)} time steps")
        
        # Perform polynomial regression to predict deorbit time
        self.predict_deorbit_time()
    
    def predict_deorbit_time(self):
        """Use polynomial regression to predict when the MHD sprint satellite will deorbit."""
        print(f"\n--- MHD Sprint Satellite Deorbit Prediction Analysis ---")
        
        if not self.time_data or len(self.time_data) < 10:
            print("Insufficient data for deorbit prediction.")
            return
        
        # Calculate altitudes over time for MHD satellite specifically
        time_hours = np.array(self.time_data) / 3600
        altitudes = np.array([np.linalg.norm(state[:3]) - 6371000 for state in self.mhd_states]) / 1000  # km
        
        # Check if satellite has already deorbited during simulation
        min_altitude = np.min(altitudes)
        current_altitude = altitudes[-1]
        
        print(f"Current MHD satellite altitude: {current_altitude:.2f} km")
        print(f"Minimum altitude reached: {min_altitude:.2f} km")
        print(f"Deorbit boundary: {TERMINATION_ALTITUDE/1000:.1f} km")
        
        # If satellite has already gone below deorbit boundary, it has deorbited
        if min_altitude <= TERMINATION_ALTITUDE/1000:
            # Find when it first crossed the deorbit boundary
            deorbit_index = np.where(altitudes <= TERMINATION_ALTITUDE/1000)[0]
            if len(deorbit_index) > 0:
                first_deorbit_time = time_hours[deorbit_index[0]]
                print(f"MHD satellite has already deorbited!")
                print(f"Deorbit occurred at: {first_deorbit_time:.2f} hours ({first_deorbit_time/24:.2f} days)")
                print(f"Altitude at deorbit: {altitudes[deorbit_index[0]]:.2f} km")
                
                # Save deorbit data
                deorbit_data = {
                    'Metric': [
                        'Deorbit Status',
                        'Deorbit Time (hours)',
                        'Deorbit Time (days)',
                        'Altitude at Deorbit (km)',
                        'Deorbit Boundary (km)',
                        'Simulation Duration (hours)'
                    ],
                    'Value': [
                        'DEORBITED',
                        f"{first_deorbit_time:.2f}",
                        f"{first_deorbit_time/24:.2f}",
                        f"{altitudes[deorbit_index[0]]:.2f}",
                        f"{TERMINATION_ALTITUDE/1000:.1f}",
                        f"{time_hours[-1]:.2f}"
                    ]
                }
                
                deorbit_df = pd.DataFrame(deorbit_data)
                deorbit_df.to_csv('deorbit_prediction.csv', index=False)
                print(f"Deorbit data saved to deorbit_prediction.csv")
                return
        
        # Calculate velocities and accelerations for trend analysis
        velocities = []
        accelerations = []
        
        for i in range(len(self.mhd_states)):
            # Velocity magnitude
            v_mag = np.linalg.norm(self.mhd_states[i][3:]) / 1000  # km/s
            velocities.append(v_mag)
            
            # Acceleration (change in velocity between steps)
            if i > 0:
                dt = (self.time_data[i] - self.time_data[i-1]) / 3600  # hours
                dv = (velocities[i] - velocities[i-1])  # km/s
                acc = dv / dt if dt > 0 else 0  # km/s²
                accelerations.append(acc)
            else:
                accelerations.append(0)
        
        # Fit polynomial to altitude data (3rd degree)
        try:
            coeffs_alt = np.polyfit(time_hours, altitudes, 3)
            poly_alt = np.poly1d(coeffs_alt)
            
            # Find when altitude reaches deorbit boundary
            def altitude_func(t):
                return poly_alt(t) - TERMINATION_ALTITUDE/1000  # deorbit boundary
            
            # Use Newton's method to find root (when altitude = 50 km)
            from scipy.optimize import fsolve
            deorbit_time_hours = fsolve(altitude_func, time_hours[-1] + 10)[0]
            
            # Check if prediction is reasonable
            if deorbit_time_hours > time_hours[-1] and deorbit_time_hours < 10000:  # reasonable range
                print(f"Polynomial regression prediction:")
                print(f"  Altitude trend: {coeffs_alt[0]:.5f}t³ + {coeffs_alt[1]:.5f}t² + {coeffs_alt[2]:.5f}t + {coeffs_alt[3]:.5f}")
                print(f"  Predicted deorbit time: {deorbit_time_hours:.2f} hours ({deorbit_time_hours/24:.2f} days)")
                print(f"  Time to deorbit: {deorbit_time_hours - time_hours[-1]:.2f} hours")
                
                # Calculate deorbit rate
                deorbit_rate = (altitudes[0] - TERMINATION_ALTITUDE/1000) / deorbit_time_hours  # km/hour
                print(f"  Average deorbit rate: {deorbit_rate:.5f} km/hour")
                
                # Save prediction to CSV
                prediction_data = {
                    'Metric': [
                        'Current Simulation Time (hours)',
                        'Predicted Deorbit Time (hours)',
                        'Time to Deorbit (hours)',
                        'Time to Deorbit (days)',
                        'Average Deorbit Rate (km/hour)',
                        'Deorbit Boundary (km)',
                        'Current Altitude (km)',
                        'Polynomial Coefficients (a*t³ + b*t² + c*t + d)',
                        'Coefficient a',
                        'Coefficient b', 
                        'Coefficient c',
                        'Coefficient d'
                    ],
                    'Value': [
                        f"{time_hours[-1]:.2f}",
                        f"{deorbit_time_hours:.2f}",
                        f"{deorbit_time_hours - time_hours[-1]:.2f}",
                        f"{(deorbit_time_hours - time_hours[-1])/24:.2f}",
                        f"{deorbit_rate:.5f}",
                        f"{TERMINATION_ALTITUDE/1000:.1f}",
                        f"{current_altitude:.2f}",
                        f"Altitude = a*t³ + b*t² + c*t + d",
                        f"{coeffs_alt[0]:.5f}",
                        f"{coeffs_alt[1]:.5f}",
                        f"{coeffs_alt[2]:.5f}",
                        f"{coeffs_alt[3]:.5f}"
                    ]
                }
                
                prediction_df = pd.DataFrame(prediction_data)
                prediction_df.to_csv('deorbit_prediction.csv', index=False)
                print(f"  Deorbit prediction saved to deorbit_prediction.csv")
                
            else:
                print(f"Polynomial regression prediction: No reasonable deorbit time found")
                print(f"  Altitude trend: {coeffs_alt[0]:.5f}t³ + {coeffs_alt[1]:.5f}t² + {coeffs_alt[2]:.5f}t + {coeffs_alt[3]:.5f}")
                print(f"  Current altitude: {current_altitude:.2f} km")
                print(f"  Altitude change rate: {coeffs_alt[2]:.5f} km/hour")
                print(f"  Deorbit boundary: {TERMINATION_ALTITUDE/1000:.1f} km")
                print(f"  Note: Satellite may not deorbit within reasonable timeframe")
                
        except Exception as e:
            print(f"Error in polynomial regression: {e}")
            print(f"Current altitude: {current_altitude:.2f} km")
            print(f"Altitude change over simulation: {altitudes[0] - altitudes[-1]:.2f} km")
            print(f"Deorbit boundary: {TERMINATION_ALTITUDE/1000:.1f} km")

def main():
    """Main function to run the complete simulation."""
    # Create simulation instance
    sim = Simulation()
    
    # Run simulation for 1 week with 3000 discrete time steps
    t_span = (0, 3600*24)
    num_steps = int(t_span[1] / 180) 
    
    # Run simulation
    sim.run_simulation(t_span, num_steps)
    
    # Generate plots
    sim.generate_plots()
    
    # Generate summary
    sim.generate_summary()
    
    print("\nSimulation complete")

if __name__ == "__main__":
    main()