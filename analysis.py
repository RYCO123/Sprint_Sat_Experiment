"""
Analysis functions for Sprint Satellite simulation.

This module contains functions for analyzing simulation results,
generating summaries, and predicting deorbit times.
"""

import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from config import (
    MASS_1U, SOLAR_POWER_PEAK, SOLAR_EFFICIENCY, ELECTRODE_DISTANCE, GENERATOR_LENGTH, MAGNETIC_DISTANCE,
    MHD_MAGNET_STRENGTH, MHD_CONDUCTOR_RESISTIVITY, MHD_CONDUCTOR_DIAMETER, MHD_CIRCUIT_RESISTANCE,
    INITIAL_ALTITUDE, INCLINATION_DEG, TERMINATION_ALTITUDE, PLASMA_TEMPERATURE, R_EARTH
)
import os


def generate_summary(time_data, mhd_states, solar_states, mhd_power, solar_power, mhd_energy, solar_energy):
    """
    Generate simulation summary and save to CSV.
    
    Args:
        time_data: List of simulation times
        mhd_states: List of MHD satellite states
        solar_states: List of solar satellite states
        mhd_power: List of MHD satellite power values
        solar_power: List of solar satellite power values
        mhd_energy: List of MHD satellite energy values
        solar_energy: List of solar satellite energy values
    """
    print("\n--- Simulation Summary ---")
    
    if not time_data:
        print("No simulation data available.")
        return
    
    # Ensure results directory exists
    results_dir = 'results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Calculate final metrics
    final_time = time_data[-1]
    mission_duration_hours = final_time / 3600
    final_mhd_state = mhd_states[-1]
    final_solar_state = solar_states[-1]
    final_mhd_alt = np.linalg.norm(final_mhd_state[:3]) - R_EARTH
    final_solar_alt = np.linalg.norm(final_solar_state[:3]) - R_EARTH
    avg_mhd_power = np.mean(mhd_power)
    avg_solar_power = np.mean(solar_power)
    total_mhd_energy = mhd_energy[-1] / 3600  # Wh
    total_solar_energy = solar_energy[-1] / 3600  # Wh

    # Section 1: General Experiment Results
    results_rows = [
        ["General Experiment Results", "", ""],
        ["Mission Duration (hours)", f"{mission_duration_hours:.5f}", f"{mission_duration_hours:.5f}"],
        ["Initial Altitude (km)", f"{INITIAL_ALTITUDE/1000:.5f}", f"{INITIAL_ALTITUDE/1000:.5f}"],
        ["Final Altitude (km)", f"{final_mhd_alt/1000:.5f}", f"{final_solar_alt/1000:.5f}"],
        ["Average Power Generated (W)", f"{avg_mhd_power:.5f}", f"{avg_solar_power:.5f}"],
        ["Total Energy Generated (Wh)", f"{total_mhd_energy:.5f}", f"{total_solar_energy:.5f}"],
        ["Power Ratio (MHD/Solar)", f"{avg_mhd_power/avg_solar_power:.5f}x", "1.00000x"],
        ["Energy Ratio (MHD/Solar)", f"{total_mhd_energy/total_solar_energy:.5f}x", "1.00000x"],
        ["", "", ""]
    ]

    # Section 2: Configuration Parameters
    config_rows = [
        ["Configuration Parameters", "", ""],
        ["Mass (kg)", f"{MASS_1U}", f"{MASS_1U}"],
        ["Electrode Distance (m)", f"{ELECTRODE_DISTANCE}", "-"],
        ["Generator Length (m)", f"{GENERATOR_LENGTH}", "-"],
        ["Magnetic Distance (m)", f"{MAGNETIC_DISTANCE}", "-"],
        ["Magnet Strength (T)", f"{MHD_MAGNET_STRENGTH}", "-"],
        ["Conductor Resistivity (Ohm*m)", f"{MHD_CONDUCTOR_RESISTIVITY}", "-"],
        ["Conductor Diameter (m)", f"{MHD_CONDUCTOR_DIAMETER}", "-"],
        ["Circuit Resistance (Ohm)", f"{MHD_CIRCUIT_RESISTANCE}", "-"],
        ["Plasma Temperature (K)", f"{PLASMA_TEMPERATURE}", "-"],
        ["Peak Power (W)", "-", f"{SOLAR_POWER_PEAK}"],
        ["Efficiency", "-", f"{SOLAR_EFFICIENCY}"],
        ["Initial Altitude (m)", f"{INITIAL_ALTITUDE}", f"{INITIAL_ALTITUDE}"],
        ["Inclination (deg)", f"{INCLINATION_DEG}", f"{INCLINATION_DEG}"],
        ["Deorbit Altitude (m)", f"{TERMINATION_ALTITUDE}", f"{TERMINATION_ALTITUDE}"],
    ]

    # Combine all rows
    all_rows = results_rows + config_rows
    summary_df = pd.DataFrame(all_rows, columns=["Metric", "MHD Sprint Satellite", "Standard Solar Satellite"])
    summary_path = os.path.join(results_dir, 'simulation_summary.csv')
    summary_df.to_csv(summary_path, index=False)
    print(summary_df.to_string(index=False))
    print(f"\nSummary saved to {summary_path}")
    
    # Print key insights
    print(f"\n--- Key Insights ---")
    print(f"MHD Sprint Satellite generates {avg_mhd_power/avg_solar_power:.1f}x more power")
    print(f"MHD Sprint Satellite generates {total_mhd_energy/total_solar_energy:.1f}x more total energy")
    print(f"Simulation completed in {mission_duration_hours:.1f} hours with {len(time_data)} time steps")


def predict_deorbit_time(time_data, mhd_states):
    """
    Use polynomial regression to predict when the MHD sprint satellite will deorbit.
    
    Args:
        time_data: List of simulation times
        mhd_states: List of MHD satellite states
    """
    print(f"\n--- MHD Sprint Satellite Deorbit Prediction Analysis ---")
    
    if not time_data or len(time_data) < 10:
        print("Insufficient data for deorbit prediction.")
        return
    
    # Calculate altitudes over time for MHD satellite specifically
    time_hours = np.array(time_data) / 3600
    altitudes = np.array([np.linalg.norm(state[:3]) - R_EARTH for state in mhd_states]) / 1000  # km
    
    min_altitude = np.min(altitudes)
    current_altitude = altitudes[-1]
    print(f"Current MHD satellite altitude: {current_altitude:.2f} km")
    print(f"Minimum altitude reached: {min_altitude:.2f} km")
    print(f"Deorbit boundary: {TERMINATION_ALTITUDE/1000:.1f} km")

    results_dir = 'results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    prediction_path = os.path.join(results_dir, 'deorbit_prediction.csv')

    # If satellite has already gone below deorbit boundary, report when
    if min_altitude <= TERMINATION_ALTITUDE/1000:
        deorbit_index = np.where(altitudes <= TERMINATION_ALTITUDE/1000)[0]
        if len(deorbit_index) > 0:
            first_deorbit_time = time_hours[deorbit_index[0]]
            print(f"MHD satellite has already deorbited!")
            print(f"Deorbit occurred at: {first_deorbit_time:.2f} hours ({first_deorbit_time/24:.2f} days)")
            print(f"Altitude at deorbit: {altitudes[deorbit_index[0]]:.2f} km")
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
            deorbit_df.to_csv(prediction_path, index=False)
            print(f"Deorbit data saved to {prediction_path}")
            return

    # Otherwise, fit a cubic polynomial and find when it will cross the boundary
    coeffs_alt = np.polyfit(time_hours, altitudes, 3)
    poly_alt = np.poly1d(coeffs_alt)
    def altitude_func(t):
        return poly_alt(t) - TERMINATION_ALTITUDE/1000

    from scipy.optimize import fsolve
    guess = max(time_hours[-1] + 10, 1.0)
    roots = fsolve(altitude_func, [guess, guess+100, guess+1000, guess+10000])
    future_roots = [r for r in roots if np.isreal(r) and r > time_hours[-1]]
    if future_roots:
        deorbit_time_hours = float(np.min(future_roots))
        print(f"Cubic regression prediction:")
        print(f"  Altitude trend: {coeffs_alt[0]:.5f}t³ + {coeffs_alt[1]:.5f}t² + {coeffs_alt[2]:.5f}t + {coeffs_alt[3]:.5f}")
        print(f"  Predicted deorbit time: {deorbit_time_hours:.2f} hours ({deorbit_time_hours/24:.2f} days)")
        print(f"  Time to deorbit: {deorbit_time_hours - time_hours[-1]:.2f} hours")
        deorbit_rate = (altitudes[0] - TERMINATION_ALTITUDE/1000) / deorbit_time_hours  # km/hour
        print(f"  Average deorbit rate: {deorbit_rate:.5f} km/hour")
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
        prediction_df.to_csv(prediction_path, index=False)
        print(f"  Deorbit prediction saved to {prediction_path}")
    else:
        print("Not enough data to find deorbit time.")
        prediction_data = {
            'Metric': ['Deorbit Prediction'],
            'Value': ['Not enough data to find deorbit time']
        }
        prediction_df = pd.DataFrame(prediction_data)
        prediction_df.to_csv(prediction_path, index=False)
        print(f"  Deorbit prediction saved to {prediction_path}") 