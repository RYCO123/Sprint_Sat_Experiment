# --- Physical Constants ---
G = 6.67430e-11  # Gravitational constant (m³/kg/s²)
M_EARTH = 5.972e24  # Earth's mass (kg)
R_EARTH = 6371000  # Earth's radius (m)
MU = G * M_EARTH  # Earth's gravitational parameter (m³/s²)

# --- 1U CubeSat Parameters (Fixed) ---
MASS_1U = 1.33  # kg
AREA_1U = 0.01  # m² (cross-sectional area)
DRAG_COEFFICIENT = 2.2

# --- MHD Sprint Satellite Parameters (Deployable Tether) ---
TETHER_LENGTH = 50.0  # m (deployable tether length - reasonably small)
TETHER_DIAMETER = 0.001  # m (1 mm diameter - reasonably small)
TETHER_RESISTIVITY = 1.68e-8  # Ω⋅m (copper resistivity)
CIRCUIT_RESISTANCE = 10.0  # Ohms (additional circuit resistance)

# --- Solar Satellite Parameters ---
SOLAR_POWER_AVERAGE = 2.5  # W (average for 1U satellite)
SOLAR_POWER_PEAK = 10.0  # W (peak power)
SOLAR_EFFICIENCY = 0.25  # 25% average efficiency

# --- Orbital Parameters ---
INITIAL_ALTITUDE = 400.0e3  # m (400 km)
INCLINATION_DEG = 20.0  # degrees
TERMINATION_ALTITUDE = 150.0e3  # m (150 km)

# --- Simulation Parameters ---
TIME_STEP = 10.0  # s
SIMULATION_START_DATE = "2024-03-01"

# --- API Caching ---
CACHE_EXPIRATION_SECONDS = 600 