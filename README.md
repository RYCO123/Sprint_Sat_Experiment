# Sprint Satellite Simulation

## Project Overview

The Sprint Satellite project demonstrates the viability of a novel 1U CubeSat concept that generates significantly more power than traditional solar-powered satellites, at the cost of a drastically reduced orbital lifetime. This "Sprint" satellite uses an electrodynamic tether (MHD - MagnetoHydroDynamic) to generate power through Earth's magnetic field, enabling high-power missions with short durations.

### Key Concept

Traditional satellites must balance power generation with orbital longevity. The Sprint Satellite concept flips this paradigm:
- **High Power**: Generates 10-100x more power than solar panels
- **Short Lifetime**: Deorbits in days/weeks instead of years
- **Perfect for**: High-power, short-duration missions like:
  - Rapid technology demonstrations
  - High-bandwidth data transmission
  - Intensive scientific measurements
  - Emergency communications

## Physics Models

### MHD "Sprint" Satellite
- **Power Generation**: Electrodynamic tether interacting with Earth's magnetic field
- **Key Equations**:
  - Motional EMF: V = (v × B) · L
  - Circuit Current: I = V/R
  - Power: P = I × V
  - Drag Force: F = I × (L × B)
- **Primary Drag**: Electrodynamic drag from power generation
- **Lifetime**: Days to weeks (depending on altitude)

### Standard Solar Satellite
- **Power Generation**: Solar panels with eclipse considerations
- **Primary Drag**: Atmospheric drag
- **Lifetime**: Years to decades
- **Power Output**: Limited by panel area and efficiency

## Features

- **Real-time API Integration**: Uses NASA CCMC and BGS IGRF APIs for accurate atmospheric and magnetic field data
- **Comprehensive Physics**: Full orbital mechanics with realistic perturbations
- **Progress Tracking**: Real-time simulation progress with live metrics
- **Data Visualization**: Three detailed plots showing power, altitude, and energy comparisons
- **Caching System**: Optimized API calls with intelligent caching

## Installation

### Prerequisites
- Python 3.7 or higher
- Internet connection (for API data)

### Setup
1. Clone the repository:
```bash
git clone https://github.com/yourusername/sprint-satellite.git
cd sprint-satellite
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Basic Simulation
Run the complete simulation:
```bash
python main.py
```

This will:
1. Simulate both satellites from 400km altitude until the MHD satellite deorbits
2. Display real-time progress with current altitudes and mission day
3. Generate three plots:
   - `altitude_decay_api.png`: Orbital decay comparison
   - `power_comparison_detail_api.png`: Power generation comparison
   - `total_energy_api.png`: Cumulative energy generation
4. Create `simulation_summary_api.csv` with key metrics

### Expected Results

The simulation demonstrates:
- **MHD Satellite**: High power (10-50W), rapid deorbit (days)
- **Solar Satellite**: Low power (1-3W), slow deorbit (years)
- **Energy Advantage**: MHD satellite generates more total energy in its short lifetime

## Configuration

Key parameters can be modified in `main.py`:

### Satellite Parameters
```python
MASS = 1.33                    # kg
AREA_CROSS_SECTION = 0.01      # m² (1U CubeSat)
DRAG_COEFFICIENT = 2.2         # dimensionless
```

### MHD Satellite
```python
TETHER_LENGTH = 100.0          # m
CIRCUIT_RESISTANCE = 50.0      # Ohms
```

### Solar Satellite
```python
SOLAR_POWER_PEAK = 10.0        # W
SOLAR_ETA_SA = 0.25           # efficiency factor
```

### Orbital Parameters
```python
INITIAL_ALTITUDE = 400.0e3     # m (400 km)
INCLINATION_DEG = 20.0         # degrees
TERMINATION_ALTITUDE = 150.0e3 # m (150 km)
```

## API Dependencies

The simulation uses three external APIs:

1. **NASA CCMC NRLMSISE-00**: Atmospheric density model
2. **BGS IGRF**: Earth's magnetic field model
3. **CelesTrak**: Solar weather indices (F10.7, Ap)

All APIs are free and publicly available. The simulation includes fallback mechanisms if APIs are unavailable.

## Output Files

### Plots
- `altitude_decay_api.png`: Shows how both satellites lose altitude over time
- `power_comparison_detail_api.png`: Compares instantaneous power generation
- `total_energy_api.png`: Shows cumulative energy generated over mission lifetime

### Data
- `simulation_summary_api.csv`: Summary table with key metrics:
  - Mission duration
  - Average power generated
  - Total energy generated
  - Initial and final altitudes

## Technical Details

### Orbital Mechanics
- Two-body gravitational model with perturbations
- Real-time atmospheric density calculation
- Magnetic field vector computation
- Eclipse detection for solar satellite

### Power Generation Models
- **MHD**: Electrodynamic tether with realistic circuit resistance
- **Solar**: Panel efficiency with sun angle and eclipse effects

### Coordinate Systems
- ECI (Earth-Centered Inertial) for orbital calculations
- NED (North-East-Down) for magnetic field data
- Proper coordinate transformations implemented

## Contributing

Contributions are welcome! Areas for improvement:
- Additional satellite configurations
- More sophisticated orbital models
- Enhanced visualization options
- Performance optimizations

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- NASA CCMC for atmospheric density data
- British Geological Survey for magnetic field data
- CelesTrak for space weather indices
- The open-source scientific Python community

## Citation

If you use this simulation in your research, please cite:
```
Sprint Satellite Simulation: Viability Study of High-Power, Short-Duration CubeSat Missions
[Your Name], [Year]
```

## Contact

For questions or contributions, please open an issue on GitHub or contact [your-email@domain.com]. 