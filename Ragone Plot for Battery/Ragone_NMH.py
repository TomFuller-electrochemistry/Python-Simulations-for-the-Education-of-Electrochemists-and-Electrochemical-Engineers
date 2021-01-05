# Ragone Plot for NMH battery
import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pyplot as plt
# Input parameters
Rint = 0.002 # internal resistance in ohm
K = 0.003 # parameter for voltage model ohm
A = -0.04 # parameter for voltage model V
B = 3     # parameter for voltage model (Ah)-1
U = 1.35  # equilibrium voltage in V
Q = 4.95  # capacity in Ah
Vco = 0.9 #cutoff voltage in V
mass = 105  #grams
# C-rates for the simulation
C = [.25,1,3,5,10,13]
#
# initialize vectors
tend = np.zeros(6)
Energy = np.zeros(6)
Edensity = np.zeros(6)
Pdensity = np.zeros(6)

# Define voltage function
def voltage(t,I):
    cap = I*t/3600
    V = U - I*Rint - K*(Q/(Q-cap))*I - A*np.exp(-cap*B)
    return V
# Function to solve for the time required to reach the cutoff voltage
def cutoff(t,I):
    F = 0.
    F = voltage(t,I) - Vco
#    print(i,t,I,F)
    return F

# Compute time at cutoff voltage for each C-rate
for i in range (0,6): 
    tguess = 0.9*1/C[i]*3600  # in seconds
    I = Q*C[i]
    tend[i] = opt.fsolve(cutoff,tguess,args=I)

# Integrate to obtain energy I*Vdt
for i in range (0,6):
    I = Q*C[i]
    ts = tend[i]/3600  # in hours
    Energy[i] = I*integrate.quad(voltage,0,ts,args=I)[0]  #  Wh
    Edensity[i] = Energy[i]/mass*1000  # kWh/kg
    Pdensity[i] = Edensity[i]/ts  #kW/kg
    
plt.plot(Pdensity,Edensity, 'ko-') # linestyle: k is black, o is symbol
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Power Density kW/kg') 
plt.ylabel('Energy Density kWh/kg')
plt.minorticks_on()
plt.grid(axis='both', which='both')
plt.show()

