# Sprint Satellite Simulation

A comprehensive orbital simulation comparing an MHD (Magnetohydrodynamic) Sprint Satellite with a standard solar-powered satellite. The MHD Sprint Satellite uses a compact MHD generator integrated into a 1U CubeSat to extract power from ionospheric plasma flow while generating orbital drag for controlled deorbiting.

---
<img width="904" height="529" alt="diagram" src="https://github.com/user-attachments/assets/65b84579-3bff-4319-abef-bb783bf96ce6" />

## System Overview

This project simulates two satellites in Low Earth Orbit (LEO):

1. **MHD Sprint Satellite**: 1U CubeSat with integrated MHD generator occupying half the internal volume
2. **Standard Solar Satellite**: Traditional solar-powered 1U CubeSat for comparison

---
 
## Key Assumptions

Deorbit Condition is reached when altitude < 60km.
Both Sprint and Standard Satellites have the same mass.
The Standard Satellite experiences no drag, the Spring Satellite only experiences MHD drag
Ionosphere conditions and solar activity were both modeled with python libraries
Time step of the ODE solver is approximately one step every 3 minutes, benchmarked against a circular orbit which maintained reasonably for a simulated month
B field of the Spring Satellite is the sum of the induced magnetic field of the Satellite + the vector of earths magnetic field that is relevant to the MHD
The system is defined for the B field to be the nob. Increasing it increases power produciton but increases drag. Models could alos be developed, changing electrode length to adjust power.

---

## Project Structure

The codebase is organized for clarity and maintainability:

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
```

---

## 1U CubeSat MHD Generator Design

The MHD Sprint Satellite is designed as a **1U CubeSat (10×10×11.35 cm)** with the internal volume divided equally:

### Internal Layout

* **50% MHD Generator**: Occupies half the internal volume (5×10×10 cm)
* **50% Payload**: Remaining half for scientific instruments, communications, and control systems

### MHD Generator Configuration

The generator is a **compact rectangular design** optimized for the 1U form factor:

* **Electrode Distance**: 0.1 m (10 cm) - distance between perpendicular electrodes
* **Generator Length**: 0.1135 m (10 cm) - longitudinal distance of the generator
* **Magnetic Distance**: 0.05 m (5 cm) - distance between magnets creating B field
* **Magnetic Field**: variable and demonstrated up to 0.01T (10,000 μT) - created with permanent or electromagnets
    * Using electromagnets allows the B field to be varied in orbit, providing active control over power generation and orbital altitude (drag) during the mission.
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

1. **Clone the repository**
   ```
   git clone https://github.com/RYCO123/Sprint_Sat_Experiment.git
   cd Sprint_Sat_Experiment
   ```

2. **Install dependencies**  
   (Recommended: use a virtual environment)
   ```
   pip install -r requirements.txt
   ```

3. **Run the simulation**
   ```
   python main.py
   ```

This will:
- Simulate both satellites for the specified duration
- Generate comparison plots and summary CSV files for each run
- Predict deorbit time for the MHD satellite

---

## Experimental Runs and Results

Below are the results for each run, with plots and summary statistics. The number in each run (e.g., 0.4x, 2.4x, etc.) indicates the ratio of MHD to solar power generation.

### 0.4x Power Run

![Comprehensive Comparison Plot](plots_0-4x/comprehensive_comparison.png)

- **Average Power Generated (W):** MHD: 0.01961, Solar: 0.04630
- **Total Energy Generated (Wh):** MHD: 0.47081, Solar: 1.11345
- **Power Ratio (MHD/Solar):** 0.42x
- **Energy Ratio (MHD/Solar):** 0.42x
- **Deorbit prediction:** ~67.2 days to deorbit from 400 km (average deorbit rate: 0.21 km/hour)
- **Power output was more stable over time compared to the solar satellite, despite lower magnitude.**

### 2.4x Power Run

![Comprehensive Comparison Plot](plots_2-4x/comprehensive_comparison.png)

- **Average Power Generated (W):** MHD: 0.10952, Solar: 0.04630
- **Total Energy Generated (Wh):** MHD: 2.62880, Solar: 1.11345
- **Power Ratio (MHD/Solar):** 2.37x
- **Energy Ratio (MHD/Solar):** 2.36x
- **Deorbit prediction:** ~36.2 days to deorbit from 400 km (average deorbit rate: 0.38 km/hour)
- **Power output was consistently stable, with moderate altitude loss due to drag.**

### 7.6x Power Run

![Comprehensive Comparison Plot](plots_7-6x/comprehensive_comparison.png)

- **Average Power Generated (W):** MHD: 0.35221, Solar: 0.04630
- **Total Energy Generated (Wh):** MHD: 8.45402, Solar: 1.11345
- **Power Ratio (MHD/Solar):** 7.61x
- **Energy Ratio (MHD/Solar):** 7.59x
- **Deorbit prediction:** ~21.8 days to deorbit from 400 km (average deorbit rate: 0.62 km/hour)
- **Power output was highly stable, with accelerated orbital decay.**

### 30.6x Power Run

![Comprehensive Comparison Plot](plots_30-6x/comprehensive_comparison.png)

- **Average Power Generated (W):** MHD: 1.41778, Solar: 0.04630
- **Total Energy Generated (Wh):** MHD: 34.03123, Solar: 1.11345
- **Power Ratio (MHD/Solar):** 30.62x
- **Energy Ratio (MHD/Solar):** 30.56x
- **Deorbit prediction:** ~9.6 days to deorbit from 400 km (average deorbit rate: 1.34 km/hour)
- **Power output was extremely stable, with the fastest deorbit rate observed.**

---

## Key Results

- **Power Generation:** The MHD Sprint Satellite can generate significantly more power than a standard solar satellite, depending on configuration. Even at lower power, its output is much more stable over time, as it is not solely dependent on direct sunlight. Instead, it relies on ionospheric plasma flow, which is less affected by orbital position relative to the sun (though still influenced by ionospheric conditions).
- **Orbital Decay:** The increased power comes at the cost of much greater drag, causing the MHD satellite to deorbit far faster than a solar satellite. The higher the power output, the faster the deorbit.
- **Stability:** Power generation for the MHD satellite is far more stable than for a solar satellite, which is subject to eclipses and sun angle. This makes the MHD approach attractive for missions requiring consistent power delivery, at the expense of mission duration.
- **Tradeoff:** There is a clear tradeoff between power generation and orbital lifetime. The MHD satellite can be tuned for more power or longer life, but not both.

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
- Polynomial regression for deorbit prediction

### Unique Design Features

- **Deorbit Time**: Predicted to be ~11.35 days (with less powerful magnets, lifetime increases while power decreases)
- **Precise Deorbit**: If using electromagnets to better control the B-field, an incredibly controlled descent can be achieved; a strong B-field creates large drag acting as an electromagnetic braking system

---

## About

An experiment to determine the viability of an MHD sprint satellite concept. A satellite that could generate significantly more power for the same form factor while generating increased drag, decreasing orbital time.

---

## License

MIT license

---

If you have any questions or want to contribute, please see [CONTRIBUTING.md](CONTRIBUTING.md).

---

This README reflects the new modular structure and preserves all technical and theoretical details from the original documentation.  
For the latest code and updates, visit the [GitHub repository](https://github.com/RYCO123/Sprint_Sat_Experiment).
