# Sprint Satellite Simulation

A comprehensive orbital simulation comparing an MHD (Magnetohydrodynamic) Sprint Satellite with a standard solar-powered satellite. The MHD Sprint Satellite uses a compact MHD generator integrated into a 1U CubeSat to extract power from ionospheric plasma flow while generating orbital drag for controlled deorbiting.

**GitHub Repository:**  
[https://github.com/RYCO123/Sprint_Sat_Experiment](https://github.com/RYCO123/Sprint_Sat_Experiment)

---

## System Overview

This project simulates two satellites in Low Earth Orbit (LEO):

1. **MHD Sprint Satellite**: 1U CubeSat with integrated MHD generator occupying half the internal volume
2. **Standard Solar Satellite**: Traditional solar-powered 1U CubeSat for comparison

---

## Project Structure

The codebase is now organized for clarity and maintainability:

```
main.py                # Entry point: runs the full simulation, analysis, and plotting
config.py              # All physical constants and configuration values
utils.py               # Utility functions (space weather, orbit, sun vector, eclipse)
physics.py             # Physics calculations (orbital, MHD, plasma, ionosphere)
satellite.py           # Satellite class (handles both MHD and solar satellites)
plotting.py            # All plotting and visualization functions
analysis.py            # Analysis, summary statistics, and deorbit prediction
simulation.py          # Core simulation orchestration logic
requirements.txt       # Python dependencies
README.md              # This documentation
LICENSE, CHANGELOG.md, CONTRIBUTING.md, .gitignore, etc.
plots/                 # Output plots
```

---

## 1U CubeSat MHD Generator Design

The MHD Sprint Satellite is designed as a **1U CubeSat (10×10×10 cm)** with the internal volume divided equally:

### Internal Layout

* **50% MHD Generator**: Occupies half the internal volume (5×10×10 cm)
* **50% Payload**: Remaining half for scientific instruments, communications, and control systems

### MHD Generator Configuration

The generator is a **compact rectangular design** optimized for the 1U form factor:

* **Electrode Distance**: 0.1 m (10 cm) - distance between perpendicular electrodes
* **Generator Length**: 0.1 m (10 cm) - longitudinal distance of the generator
* **Magnetic Distance**: 0.05 m (5 cm) - distance between magnets creating B field
* **Magnetic Field**: 0.01 Tesla (10,000 μT) - enhanced field strength from permanent magnets
* **Electrodes**: Conductive electrodes for current collection
* **Orientation**: Electrodes automatically align perpendicular to velocity and magnetic field

### Physics Implementation

The MHD generator operates based on fundamental electromechanical principles from MIT's Electromechanical Dynamics:

#### Motional EMF (Electromotive Force)

From Woodson & Melcher's treatment of moving media (Chapter 6), the motional EMF is:

```
V_emf = (v × B) · L
```

Where:

* `v` = satellite velocity vector (m/s)
* `B` = total magnetic field vector (Tesla)
* `L` = electrode distance vector (m)

#### Current Flow and Power Generation

Following the electromechanical coupling principles (Chapter 3):

```
I = V_emf / R_total
P_generated = V_emf × I
```

Where `R_total` includes:

* Conductor resistance (copper wiring)
* Plasma resistance (ionospheric conductivity)
* Circuit resistance (additional losses)

#### Lorentz Drag Force

The braking force causing orbital decay follows the Lorentz force density (Chapter 8):

```
F_drag = I × (L × B)
```

### Plasma Properties

* **Density**: International Reference Ionosphere (IRI) model via PyIRI library
* **Conductivity**: Spitzer conductivity formula with electron-ion collisions
* **Temperature**: 1000K typical ionospheric temperature
* **Magnetic Field**: PyIRI magnetic field parameters (inclination, dip angle, magnetic latitude)

---

## Technical Specifications

### Orbital Parameters

* **Initial Altitude**: 400 km
* **Inclination**: 20°
* **Satellite Mass**: 1.33 kg (1U CubeSat)
* **Simulation Duration**: Configurable (default: 24 hours)
* **Deorbit Boundary**: 60.1 km altitude

### MHD Generator Parameters

* **Electrode Distance**: 0.1 m (10 cm)
* **Generator Length**: 0.1 m (10 cm)
* **Magnetic Distance**: 0.05 m (5 cm)
* **Magnetic Field**: 0.01 Tesla (Earth's field + artificial magnets)
* **Conductor Material**: Copper (1.68e-8 Ω⋅m resistivity)
* **Conductor Diameter**: 1 mm
* **Circuit Resistance**: 5.0 Ω

### Environmental Models

* **Magnetic Field**: PyIRI magnetic field parameters with dipole fallback
* **Plasma Density**: Direct PyIRI IRI_density_1day calls
* **Solar Activity**: F10.7 solar flux integration

---

## Usage

1. **Install dependencies**  
   (Recommended: use a virtual environment)
   ```
   pip install -r requirements.txt
   ```

2. **Run the simulation**
   ```
   python main.py
   ```

This will:

- Simulate both satellites for the specified duration
- Generate comparison plots in the `plots/` directory
- Create a summary CSV file with key metrics
- Predict deorbit time for the MHD satellite

---

## Output Files

- `simulation_summary.csv`: Key metrics comparison
- `deorbit_prediction.csv`: MHD satellite deorbit analysis
- `plots/orbital_trajectories_2d.png`: 2D orbital comparison
- `plots/orbital_trajectories_3d.png`: 3D orbital comparison
- `plots/comprehensive_comparison.png`: Power and energy analysis

---

## Key Results

The simulation demonstrates:

- **Power Generation**: MHD satellite generates significantly more power than solar (typically 22–27x, depending on parameters)
- **Orbital Decay**: MHD drag causes faster altitude loss compared to solar satellite
- **Energy Efficiency**: Higher total energy generation despite shorter mission duration
- **Controlled Deorbiting**: Predictable orbital decay for end-of-life management

---

## Theoretical Foundation

This simulation is based on the comprehensive treatment of electromechanical systems presented in:

**Woodson, Herbert H., and James R. Melcher. Electromechanical Dynamics. Part III: Elastic and Fluid Media.** Massachusetts Institute of Technology: MIT OpenCourseWare.

- <https://ocw.mit.edu/ans7870/resources/woodson/textbook/emd_part1.pdf>
- <https://ocw.mit.edu/ans7870/resources/woodson/textbook/emd_part2.pdf>
- <https://ocw.mit.edu/ans7870/resources/woodson/textbook/emd_part3.pdf>

Key theoretical concepts implemented:

1. **Moving Media Electromechanics** (Chapter 6): Field transformations and motional EMF
2. **Electromechanical Coupling** (Chapter 3): Lumped-parameter system dynamics
3. **Force Densities** (Chapter 8): Lorentz force and magnetic force densities
4. **MHD Interactions** (Chapter 13): Magnetohydrodynamic flow and power generation

The MHD generator design follows the principles of traveling-wave MHD interactions and electromechanical energy conversion as described in the MIT textbook.

---

## Technical Details

### Implementation Notes

- Modular Python package structure for clarity and maintainability
- Direct PyIRI calls for each evaluation step (no caching)
- Real-time ionospheric data integration
- Polynomial regression for deorbit prediction
- Comprehensive error handling and fallback models

### Performance Characteristics for the 0.01 Tesla magnetic field

- **MHD Power Generation**: 22–27x higher than solar satellite
- **Energy Efficiency**: 21–27x more total energy
- **Orbital Decay Rate**: ~0.89–1.25 km/hour average
- **Simulation Speed**: Real-time PyIRI integration

### Unique Design Features

- **Deorbit Time**: Predicted to be ~11.35 days (with less powerful magnets, lifetime increases while power decreases)
- **Precise Deorbit**: If using electromagnets to better control the B-field, an incredibly controlled descent can be achieved; a strong B-field creates large drag acting as an electromagnetic braking system

---

## About

An experiment to determine the viability of an MHD sprint satellite concept. A satellite that could generate significantly more power for the same form factor while generating increased drag, decreasing orbital time.

**GitHub Repository:**  
[https://github.com/RYCO123/Sprint_Sat_Experiment](https://github.com/RYCO123/Sprint_Sat_Experiment)

---

## License

MIT license

---

If you have any questions or want to contribute, please see [CONTRIBUTING.md](CONTRIBUTING.md).

---

This README reflects the new modular structure and preserves all technical and theoretical details from the original documentation.  
For the latest code and updates, visit the [GitHub repository](https://github.com/RYCO123/Sprint_Sat_Experiment).