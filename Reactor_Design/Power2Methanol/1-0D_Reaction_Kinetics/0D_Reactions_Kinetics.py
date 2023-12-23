# Import Libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

#Define the reactions system
def reaction_system(concentrations, time):

    c_CH4O, c_CO, c_CO2, c_H2, c_H2O = concentrations

    T = 550 
    
    ###Physical Constant
    R_const = 8.3144598 #Universal Gas Constant (molar gas constant) [J/(mol*K)]

    ###Partial Pressure
    P_CO = c_CO*R_const*T #CO, Partial Pressure
    P_H2 = c_H2*R_const*T #H2, Partial Pressure
    P_CH4O = c_CH4O*R_const*T #CH4O, Partial Pressure
    P_CO2 = c_CO2*R_const*T #CO2, Partial Pressure
    P_H2O =  c_H2O*R_const*T #H2O, Partial Pressur
    
    ###Rate Constants

    k1 = 4.0638e-6*np.exp(-11695/(R_const*T)) 
    k2 = 1.5188e-33*np.exp(-266010/(R_const*T)) 
    k3 = 9.0421e8*np.exp(-112860/(R_const*T)) 
    KCO = 8.3965e-11*np.exp(118270/(R_const*T)) 
    KCO2 = 1.7214e-10*np.exp(81287/(R_const*T)) 
    K0 = 3.5408e12*np.exp(19832/(R_const*T)) 
    K1 = 2.5813e10*np.exp(26788/(R_const*T)) 
    K3 =  6.1221e-13*np.exp(125226/(R_const*T)) 
    K4 =  4.3676e-12*np.exp(115080/(R_const*T)) 
    r0 = (1+KCO*P_CO+KCO2*P_CO2)*(P_H2**0.5+K4*P_H2O) 
    Keq1 = KCO/K0 
    Keq2 = KCO2/K1 
    Keq3= KCO2/K3 
    
    ###Define the rate equations
    r1 = (1/r0)*(k1*KCO*(P_CO*P_H2**1.5-P_CH4O/(P_H2**0.5*Keq1))) # [mol/gcat*s]
    r2 = (1/r0)*(k2*KCO2*(P_CO2*P_H2**1.5-P_CH4O*P_H2O/(P_H2**1.5*Keq2))) # [mol/gcat*s]
    r3 = (1/r0)*(k3*KCO2*(P_CO2*P_H2-P_CO*P_H2O/Keq3))


    dc_CH4O_dt = r1 + r2
    dc_CO2_dt = -r2 - r3
    dc_CO_dt = -r1 + r3
    dc_H2_dt = -2*r1-3*r2-r3
    dc_H2O_dt = r2 + r3

    return [ dc_CH4O_dt, dc_CO2_dt, dc_CO_dt, dc_H2_dt, dc_H2O_dt   ]




###Set up Initial Conditions
initial_concentrations = [ 0, 0, 1, 6, 0] # [mol/m^3]
#Temperature = 550 #[K]
time_points = np.linspace(0, 50, 1000) # [s]

#Solve the ODEs
concentrations = odeint(reaction_system, initial_concentrations, time_points, atol=1e-8, rtol=1e-6)

#Plot the results
plt.figure(figsize=(10, 6))

plt.plot(time_points, concentrations[:, 0], label='CH4O')
plt.plot(time_points, concentrations[:, 1], label='CO')
plt.plot(time_points, concentrations[:, 2], label='CO2')
plt.plot(time_points, concentrations[:, 3], label='H2')
plt.plot(time_points, concentrations[:, 4], label='H2O')


plt.xlabel('Time')
plt.ylabel('Concentration [mol/m^3]')
plt.legend()
plt.title('Reaction System')

plt.show()
