import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from constants import R_EARTH

def plot_2D_orbit_animation(x_orbit, y_orbit, satellite_height, R_EARTH):
    """
    Plots and animates a satellite's orbit.

    Args:
        x_orbit (array-like): X coordinates of the orbit.
        y_orbit (array-like): Y coordinates of the orbit.
        satellite_height (float): Initial orbital height (in meters).
        R_EARTH (float): Radius of Earth (in meters).
    """
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_title('Satellite Orbit Animation')
    ax.set_xlabel('X Position (meters)')
    ax.set_ylabel('Y Position (meters)')
    ax.grid(True, linestyle='--', alpha=0.6)

    earth_circle = plt.Circle((0, 0), R_EARTH, color='deepskyblue', zorder=5)
    ax.add_patch(earth_circle)

    ax.plot(x_orbit, y_orbit, 'w--', alpha=0.8, label='Orbit Path')
    satellite_point, = ax.plot([], [], 'ro', markersize=6, label='Satellite')

    max_range = (R_EARTH + satellite_height) * 1.1
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)

    fig.set_facecolor('#1c1c1c')
    ax.set_facecolor('#2b2b2b')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.title.set_color('white')

    def init():
        satellite_point.set_data([], [])
        return satellite_point,

    def update(frame):
        satellite_point.set_data([x_orbit[frame]], [y_orbit[frame]])
        return satellite_point,

    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(x_orbit),
        init_func=init,
        blit=True,
        interval=30
    )

    ax.legend()
    plt.show()

def plot_power_comparison(time_data, mhd_power, solar_power, save_path='power_comparison.png'):
    """
    Plot power generation comparison between MHD and solar satellites.
    
    Args:
        time_data: Time array in seconds
        mhd_power: MHD satellite power array
        solar_power: Solar satellite power array
        save_path: Path to save the plot
    """
    time_hours = np.array(time_data) / 3600
    
    plt.figure(figsize=(12, 6))
    plt.plot(time_hours, mhd_power, label='MHD Sprint Satellite', linewidth=2, color='red')
    plt.plot(time_hours, solar_power, label='Standard Solar Satellite', linewidth=2, color='blue')
    
    plt.xlabel('Time (hours)')
    plt.ylabel('Power Generated (W)')
    plt.title('Power Generation Comparison: MHD vs Standard Solar')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add some statistics
    avg_mhd = np.mean(mhd_power)
    avg_solar = np.mean(solar_power)
    plt.text(0.02, 0.98, f'Avg MHD Power: {avg_mhd:.2f} W\nAvg Solar Power: {avg_solar:.2f} W\nRatio: {avg_mhd/avg_solar:.1f}x', 
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Power comparison plot saved to {save_path}")

def plot_total_energy_comparison(time_data, mhd_energy, solar_energy, save_path='total_energy_comparison.png'):
    """
    Plot total energy generation comparison between MHD and solar satellites.
    
    Args:
        time_data: Time array in seconds
        mhd_energy: MHD satellite cumulative energy array (Joules)
        solar_energy: Solar satellite cumulative energy array (Joules)
        save_path: Path to save the plot
    """
    time_hours = np.array(time_data) / 3600
    mhd_energy_wh = np.array(mhd_energy) / 3600  # Convert to Watt-hours
    solar_energy_wh = np.array(solar_energy) / 3600  # Convert to Watt-hours
    
    plt.figure(figsize=(12, 6))
    plt.plot(time_hours, mhd_energy_wh, label='MHD Sprint Satellite', linewidth=2, color='red')
    plt.plot(time_hours, solar_energy_wh, label='Standard Solar Satellite', linewidth=2, color='blue')
    
    plt.xlabel('Time (hours)')
    plt.ylabel('Total Energy Generated (Wh)')
    plt.title('Cumulative Energy Generation: MHD vs Standard Solar')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add some statistics
    final_mhd = mhd_energy_wh[-1]
    final_solar = solar_energy_wh[-1]
    plt.text(0.02, 0.98, f'Total MHD Energy: {final_mhd:.1f} Wh\nTotal Solar Energy: {final_solar:.1f} Wh\nRatio: {final_mhd/final_solar:.1f}x', 
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Total energy comparison plot saved to {save_path}")

def plot_orbital_trajectories_2d(mhd_states, solar_states, save_path='orbital_trajectories_2d.png'):
    """
    Plot 2D orbital trajectories of both satellites.
    
    Args:
        mhd_states: Array of MHD satellite state vectors
        solar_states: Array of solar satellite state vectors
        save_path: Path to save the plot
    """
    # Extract position data
    mhd_positions = np.array([state[:3] for state in mhd_states])
    solar_positions = np.array([state[:3] for state in solar_states])
    
    # Create 2D plot (X-Y plane)
    plt.figure(figsize=(12, 10))
    
    # Plot trajectories
    plt.plot(mhd_positions[:, 0]/1000, mhd_positions[:, 1]/1000, 
            label='MHD Sprint Satellite', linewidth=2, color='red')
    plt.plot(solar_positions[:, 0]/1000, solar_positions[:, 1]/1000, 
            label='Standard Solar Satellite', linewidth=2, color='blue')
    
    # Add Earth
    earth_circle = plt.Circle((0, 0), R_EARTH/1000, color='lightblue', alpha=0.7, label='Earth')
    plt.gca().add_patch(earth_circle)
    
    # Mark start and end points
    plt.plot(mhd_positions[0, 0]/1000, mhd_positions[0, 1]/1000, 'ro', markersize=8, label='MHD Start')
    plt.plot(mhd_positions[-1, 0]/1000, mhd_positions[-1, 1]/1000, 'r*', markersize=10, label='MHD End')
    plt.plot(solar_positions[0, 0]/1000, solar_positions[0, 1]/1000, 'bo', markersize=8, label='Solar Start')
    plt.plot(solar_positions[-1, 0]/1000, solar_positions[-1, 1]/1000, 'b*', markersize=10, label='Solar End')
    
    plt.xlabel('X Position (km)')
    plt.ylabel('Y Position (km)')
    plt.title('2D Orbital Trajectories (X-Y Plane)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"2D orbital trajectories plot saved to {save_path}")

def plot_orbital_trajectories_3d(mhd_states, solar_states, save_path='orbital_trajectories_3d.png'):
    """
    Plot 3D orbital trajectories of both satellites.
    
    Args:
        mhd_states: Array of MHD satellite state vectors
        solar_states: Array of solar satellite state vectors
        save_path: Path to save the plot
    """
    # Extract position data
    mhd_positions = np.array([state[:3] for state in mhd_states])
    solar_positions = np.array([state[:3] for state in solar_states])
    
    # Create 3D plot
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot trajectories
    ax.plot(mhd_positions[:, 0]/1000, mhd_positions[:, 1]/1000, mhd_positions[:, 2]/1000, 
           label='MHD Sprint Satellite', linewidth=2, color='red')
    ax.plot(solar_positions[:, 0]/1000, solar_positions[:, 1]/1000, solar_positions[:, 2]/1000, 
           label='Standard Solar Satellite', linewidth=2, color='blue')
    
    # Add Earth as a sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = R_EARTH/1000 * np.outer(np.cos(u), np.sin(v))
    y = R_EARTH/1000 * np.outer(np.sin(u), np.sin(v))
    z = R_EARTH/1000 * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color='lightblue', alpha=0.3)
    
    # Mark start and end points
    ax.scatter(mhd_positions[0, 0]/1000, mhd_positions[0, 1]/1000, mhd_positions[0, 2]/1000, 
              color='red', s=100, marker='o', label='MHD Start')
    ax.scatter(mhd_positions[-1, 0]/1000, mhd_positions[-1, 1]/1000, mhd_positions[-1, 2]/1000, 
              color='red', s=150, marker='*', label='MHD End')
    ax.scatter(solar_positions[0, 0]/1000, solar_positions[0, 1]/1000, solar_positions[0, 2]/1000, 
              color='blue', s=100, marker='o', label='Solar Start')
    ax.scatter(solar_positions[-1, 0]/1000, solar_positions[-1, 1]/1000, solar_positions[-1, 2]/1000, 
              color='blue', s=150, marker='*', label='Solar End')
    
    ax.set_xlabel('X Position (km)')
    ax.set_ylabel('Y Position (km)')
    ax.set_zlabel('Z Position (km)')
    ax.set_title('3D Orbital Trajectories')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"3D orbital trajectories plot saved to {save_path}")

def plot_altitude_comparison(time_data, mhd_states, solar_states, save_path='altitude_comparison.png'):
    """
    Plot altitude comparison between both satellites.
    
    Args:
        time_data: Time array in seconds
        mhd_states: Array of MHD satellite state vectors
        solar_states: Array of solar satellite state vectors
        save_path: Path to save the plot
    """
    time_hours = np.array(time_data) / 3600
    
    # Calculate altitudes
    mhd_altitudes = np.array([np.linalg.norm(state[:3]) - R_EARTH for state in mhd_states]) / 1000
    solar_altitudes = np.array([np.linalg.norm(state[:3]) - R_EARTH for state in solar_states]) / 1000
    
    plt.figure(figsize=(12, 6))
    plt.plot(time_hours, mhd_altitudes, label='MHD Sprint Satellite', linewidth=2, color='red')
    plt.plot(time_hours, solar_altitudes, label='Standard Solar Satellite', linewidth=2, color='blue')
    
    # Add termination altitude line
    plt.axhline(y=150, color='gray', linestyle='--', alpha=0.7, label='Termination Altitude (150 km)')
    
    plt.xlabel('Time (hours)')
    plt.ylabel('Altitude (km)')
    plt.title('Altitude Decay Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add some statistics
    initial_alt = mhd_altitudes[0]
    final_mhd_alt = mhd_altitudes[-1]
    final_solar_alt = solar_altitudes[-1]
    plt.text(0.02, 0.98, f'Initial Altitude: {initial_alt:.0f} km\nMHD Final: {final_mhd_alt:.0f} km\nSolar Final: {final_solar_alt:.0f} km', 
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Altitude comparison plot saved to {save_path}")

def create_comprehensive_comparison(time_data, mhd_states, solar_states, mhd_power, solar_power, 
                                  mhd_energy, solar_energy, save_path='comprehensive_comparison.png'):
    """
    Create a comprehensive comparison plot with all metrics.
    
    Args:
        time_data: Time array in seconds
        mhd_states: Array of MHD satellite state vectors
        solar_states: Array of solar satellite state vectors
        mhd_power: MHD satellite power array
        solar_power: Solar satellite power array
        mhd_energy: MHD satellite cumulative energy array
        solar_energy: Solar satellite cumulative energy array
        save_path: Path to save the plot
    """
    time_hours = np.array(time_data) / 3600
    mhd_energy_wh = np.array(mhd_energy) / 3600
    solar_energy_wh = np.array(solar_energy) / 3600
    
    # Calculate altitudes
    mhd_altitudes = np.array([np.linalg.norm(state[:3]) - R_EARTH for state in mhd_states]) / 1000
    solar_altitudes = np.array([np.linalg.norm(state[:3]) - R_EARTH for state in solar_states]) / 1000
    
    # Create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Altitude comparison
    ax1.plot(time_hours, mhd_altitudes, label='MHD Sprint', linewidth=2, color='red')
    ax1.plot(time_hours, solar_altitudes, label='Standard Solar', linewidth=2, color='blue')
    ax1.axhline(y=150, color='gray', linestyle='--', alpha=0.7, label='Termination')
    ax1.set_xlabel('Time (hours)')
    ax1.set_ylabel('Altitude (km)')
    ax1.set_title('Altitude Decay')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Power comparison
    ax2.plot(time_hours, mhd_power, label='MHD Sprint', linewidth=2, color='red')
    ax2.plot(time_hours, solar_power, label='Standard Solar', linewidth=2, color='blue')
    ax2.set_xlabel('Time (hours)')
    ax2.set_ylabel('Power (W)')
    ax2.set_title('Power Generation')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Energy comparison
    ax3.plot(time_hours, mhd_energy_wh, label='MHD Sprint', linewidth=2, color='red')
    ax3.plot(time_hours, solar_energy_wh, label='Standard Solar', linewidth=2, color='blue')
    ax3.set_xlabel('Time (hours)')
    ax3.set_ylabel('Energy (Wh)')
    ax3.set_title('Cumulative Energy')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: 2D trajectory
    mhd_positions = np.array([state[:3] for state in mhd_states])
    solar_positions = np.array([state[:3] for state in solar_states])
    ax4.plot(mhd_positions[:, 0]/1000, mhd_positions[:, 1]/1000, 
            label='MHD Sprint', linewidth=2, color='red')
    ax4.plot(solar_positions[:, 0]/1000, solar_positions[:, 1]/1000, 
            label='Standard Solar', linewidth=2, color='blue')
    earth_circle = plt.Circle((0, 0), R_EARTH/1000, color='lightblue', alpha=0.7)
    ax4.add_patch(earth_circle)
    ax4.set_xlabel('X Position (km)')
    ax4.set_ylabel('Y Position (km)')
    ax4.set_title('2D Trajectories')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.axis('equal')
    
    plt.suptitle('Sprint Satellite vs Standard Solar Satellite Comparison', fontsize=16)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Comprehensive comparison plot saved to {save_path}")
