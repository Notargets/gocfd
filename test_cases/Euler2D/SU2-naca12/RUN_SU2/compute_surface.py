import pandas as pd
import numpy as np

# Read the CSV file
df = pd.read_csv('surface_flow.csv')

# Freestream conditions
gamma = 1.4
M_inf = 0.8
alpha = 1.25 * np.pi / 180.0

# SU2 uses REF_DIMENSIONALIZATION= FREESTREAM_PRESS_EQ_ONE
# This means: p_inf = 1.0 (not 1/gamma)
p_inf = 1.0

# For this non-dimensionalization with p_inf = 1:
# From ideal gas: p = rho * R * T, and a^2 = gamma * R * T
# So: p_inf = rho_inf * R * T_inf = 1
# And: a_inf^2 = gamma * R * T_inf = gamma * (p_inf / rho_inf) = gamma / rho_inf
# Since V_inf = M_inf * a_inf, we have: V_inf^2 = M_inf^2 * gamma / rho_inf

# For REF_DIMENSIONALIZATION= FREESTREAM_PRESS_EQ_ONE:
rho_inf = gamma  # This makes p_inf = 1 and gives correct relations
a_inf = 1.0      # Speed of sound
V_inf = M_inf * a_inf  # Freestream velocity magnitude

# Freestream dynamic pressure
q_inf = 0.5 * rho_inf * V_inf**2

print(f"SU2 Non-dimensionalization (FREESTREAM_PRESS_EQ_ONE):")
print(f"rho_inf: {rho_inf:.6f}")
print(f"p_inf: {p_inf:.6f}")
print(f"V_inf: {V_inf:.6f}")
print(f"q_inf: {q_inf:.6f}")

# Extract variables
x = df['x']
rho = df['Density']
rho_u = df['Momentum_x']
rho_v = df['Momentum_y']
rho_E = df['Energy']

# Calculate primitive variables
u = rho_u / rho
v = rho_v / rho
V_mag = np.sqrt(u**2 + v**2)

# Calculate pressure
p = (gamma - 1) * (rho_E - 0.5 * rho * V_mag**2)

# Calculate Mach number
a = np.sqrt(gamma * p / rho)
Mach = V_mag / a

# Calculate Cp using SU2's non-dimensionalization
Cp = (p - p_inf) / q_inf

print(f"Pressure range: {p.min():.6f} to {p.max():.6f}")
print(f"Cp range: {Cp.min():.3f} to {Cp.max():.3f}")
print(f"Mach range: {Mach.min():.3f} to {Mach.max():.3f}")

# Create output dataframe
output = pd.DataFrame({
    'x': x,
    'Mach': Mach,
    'Cp': Cp
})

# Save to file
output.to_csv('mach_cp_data.dat', index=False, header=False)
print("Data saved to mach_cp_data.dat")
