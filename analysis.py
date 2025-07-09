"""
Analysis functions for Sprint Satellite simulation.

This module contains functions for analyzing simulation results,
generating summaries, and predicting deorbit times.
"""

import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from config import R_EARTH, TERMINATION_ALTITUDE, INITIAL_ALTITUDE


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
    
    # Calculate final metrics
    final_time = time_data[-1]
    mission_duration_hours = final_time / 3600
    
    # Final altitudes
    final_mhd_state = mhd_states[-1]
    final_solar_state = solar_states[-1]
    final_mhd_alt = np.linalg.norm(final_mhd_state[:3]) - R_EARTH
    final_solar_alt = np.linalg.norm(final_solar_state[:3]) - R_EARTH
    
    # Average powers
    avg_mhd_power = np.mean(mhd_power)
    avg_solar_power = np.mean(solar_power)
    
    # Total energies
    total_mhd_energy = mhd_energy[-1] / 3600  # Convert to Wh
    total_solar_energy = solar_energy[-1] / 3600  # Convert to Wh
    
    # Create summary dataframe
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
            f"{INITIAL_ALTITUDE/1000:.5f}",
            f"{final_mhd_alt/1000:.5f}",
            f"{avg_mhd_power:.5f}",
            f"{total_mhd_energy:.5f}",
            f"{avg_mhd_power/avg_solar_power:.5f}x",
            f"{total_mhd_energy/total_solar_energy:.5f}x"
        ],
        'Standard Solar Satellite': [
            f"{mission_duration_hours:.5f}",
            f"{INITIAL_ALTITUDE/1000:.5f}",
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
    
    for i in range(len(mhd_states)):
        # Velocity magnitude
        v_mag = np.linalg.norm(mhd_states[i][3:]) / 1000  # km/s
        velocities.append(v_mag)
        
        # Acceleration (change in velocity between steps)
        if i > 0:
            dt = (time_data[i] - time_data[i-1]) / 3600  # hours
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