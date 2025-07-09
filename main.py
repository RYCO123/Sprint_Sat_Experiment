"""
Main entry point for Sprint Satellite simulation.

This script orchestrates the complete simulation including orbital dynamics,
power generation, analysis, and visualization.
"""

from simulation import Simulation


def main():
    """Main function to run the complete simulation."""
    # Create and run simulation
    sim = Simulation()
    t_span = (0, 3600*4) 
    num_steps = int(t_span[1] / 180) 
    
    # Run the simulation
    sim.run_simulation(t_span, num_steps)
    
    # Generate plots and analysis
    sim.generate_plots()
    sim.generate_summary()
    
    print("\nSimulation complete")


if __name__ == "__main__":
    main() 