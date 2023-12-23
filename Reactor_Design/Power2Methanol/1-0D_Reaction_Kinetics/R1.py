# Import Libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

#Define the reactions system
def reaction_system(concentrations, time):

    c_CH4O, c_CO, c_H2 = concentrations

    T = 550 #[K]
    
    ###Physical Constant
    R_const = 8.3144598 #Universal Gas Constant (molar gas constant) [J/(mol*K)]

    ###Partial Pressure
    P_CO = c_CO*R_const*T #CO, Partial Pressure
    P_H2 = c_H2*R_const*T #H2, Partial Pressure
    P_CH4O = c_CH4O*R_const*T #CH4O, Partial Pressure
    
    
    ###Rate Constants

    k1 = 4.0638e-6*np.exp(-11695/(R_const*T))
    KCO = 8.3965e-11*np.exp(118270/(R_const*T)) 
    K0 = 3.5408e12*np.exp(19832/(R_const*T))  
    Keq1 = KCO/K0 
    
    
    ###Define the rate equations
    r1 = (k1*KCO*(P_CO*P_H2**1.5-P_CH4O/(P_H2**0.5*Keq1))) # [mol/gcat*s]
    

    dc_CH4O_dt = r1 
    dc_CO_dt = -r1 
    dc_H2_dt = -2*r1
    

    return [ dc_CH4O_dt, dc_CO_dt, dc_H2_dt ]




###Set up Initial Conditions
initial_concentrations = [ 0, 1, 6] # [mol/m^3]
#Temperature = 550 #[K]
time_points = np.linspace(0, 50, 1000) # [s]

#Solve the ODEs
concentrations = odeint(reaction_system, initial_concentrations, time_points)

#Plot the results
plt.figure(figsize=(10, 6))

plt.plot(time_points, concentrations[:, 0], label='CH4O')
plt.plot(time_points, concentrations[:, 1], label='CO')
plt.plot(time_points, concentrations[:, 2], label='H2')



plt.xlabel('Time')
plt.ylabel('Concentration [mol/m^3]')
plt.legend()
plt.title('Reaction System')

plt.show()
