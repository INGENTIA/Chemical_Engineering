# Import Libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

#Define the reactions system
def reaction_system(concentrations, time):

    c_CH4O, c_CO, c_CO2, c_H2, c_H2O = concentrations

    T = 550 #[K]
    
    ###Physical Constant
    R_const = 8.3144598 #Universal Gas Constant (molar gas constant) [J/(mol*K)]

    ###Partial Pressure
    P_CO = c_CO*R_const*T #CO, Partial Pressure
    P_H2 = c_H2*R_const*T #H2, Partial Pressure
    P_CH4O = c_CH4O*R_const*T #CH4O, Partial Pressure
    P_CO2 = c_CO2*R_const*T #CO2, Partial Pressure
    P_H2O =  c_H2O*R_const*T #H2O, Partial Pressure
    
    ###Rate Constants

    kA = 5.95e9*np.exp(-113711/(R_const*T)) 
    kB = 1.16e11*np.exp(-126573/(R_const*T)) 
    kC = 2.24e6*np.exp(-68252/(R_const*T)) 
    KH2O = 3.80e-10*np.exp(80876/(R_const*T)) 
    KCO = 8.00e-6*np.exp(58015/(R_const*T)) 
    KCO2 = 1.02e-6*np.exp(67439/(R_const*T)) 
    KH2 = 2.71e2*np.exp(-6291/(R_const*T)) 
    KPA = 10**((4817/T)-9.83) 
    KPB = 10**((-2090/T)-2.018) 
    KPC = KPA*KPB 
    r0 = (1+KCO+P_CO)*(1+KH2**0.5*P_H2**0.5+KH2O*P_H2O) 
    
    ###Define the rate equations
    r1 = (1/r0)*(kA*KCO*(P_CO*P_H2**1.5-P_CH4O/(KPA*P_H2**0.5))) # [mol/gcat*s]
    r2 = (1/r0)*(kC*KCO2*(P_CO2*P_H2**1.5-(P_H2O*P_CH4O/(KPC*P_H2**1.5)))) # [mol/gcat*s]
    r3 = (1/r0)*(kB*KCO2*(P_CO2*P_H2-(P_CO*P_H2O/KPB)))


    dc_CH4O_dt = r1 + r2
    dc_CO2_dt = -r2 - r3
    dc_CO_dt = -r1 + r3
    dc_H2_dt = -2*r1-3*r2-r3
    dc_H2O_dt = r2 + r3

    return [ dc_CH4O_dt, dc_CO2_dt, dc_CO_dt, dc_H2_dt, dc_H2O_dt   ]




###Set up Initial Conditions
initial_concentrations = [ 0, 0, 1, 6, 0] # [mol/m^3]
#Temperature = 550 #[K]
time_points = np.linspace(0, 0.1, 1000) # [s]

#Solve the ODEs
concentrations = odeint(reaction_system, initial_concentrations, time_points)

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
