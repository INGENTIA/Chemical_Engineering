import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the reaction system
def reaction_system(concentrations, time):
    A, B, C = concentrations

    k1 = 0.1
    k2 = 0.02
    k3 = 0.05

    # Define the rate equations
    dA_dt = -k1 * A - 2 * k3 * A**2 * B
    dB_dt = k1 * A - k2 * B * C
    dC_dt = k2 * B * C + 3 * k3 * A**2 * B

    return [dA_dt, dB_dt, dC_dt]

# Set up the initial conditions
initial_concentrations = [1, 0, 0]
time_points = np.linspace(0, 50, 1000)

# Solve the ODEs
concentrations = odeint(reaction_system, initial_concentrations, time_points)

# Plot the results
plt.figure(figsize=(10, 6))

plt.plot(time_points, concentrations[:, 0], label='A')
plt.plot(time_points, concentrations[:, 1], label='B')
plt.plot(time_points, concentrations[:, 2], label='C')

plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()
plt.title('Reaction System')

plt.show()
