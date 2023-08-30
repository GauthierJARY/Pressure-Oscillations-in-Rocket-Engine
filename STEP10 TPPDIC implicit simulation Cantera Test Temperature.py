# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 17:26:29 2023

@author: gauth
"""

# =============================================================================
# Imported modules 
# =============================================================================

import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from scipy.optimize import fsolve
import time 
import cantera as ct 




# =============================================================================
# Geometric parametrisation 
# =============================================================================

l_tank = 10
d_tank = 1
A_tank = (d_tank/2)**2 * 3.14
V_tank = l_tank * A_tank

l_pipe = 2
d_pipe = 0.01
A_pipe = (d_pipe/2)**2 * 3.14
V_pipe = l_pipe * A_pipe

l_dome = 0.1
d_dome = 0.15
A_dome = (d_dome/2)**2 * 3.14
V_dome = l_dome * A_dome

NI = 5 # number of injectors at the end of injection dome of each species

l_injector = 0.05
d_injector = 0.003
A_injector = (d_injector/2)**2 * 3.14
V_injector = l_injector * A_injector

l_combustor = 0.2
d_combustor = 0.08 
A_combustor = (d_combustor/2)**2 * 3.14
V_combustor =  A_combustor * l_combustor 
A_throat = 0.01

# =============================================================================
# Constant values and parameters
# =============================================================================

R = 8.734
R_cc = R/(29e-3)

T3 = 3500 #K

t = 0
dt = 1e-3

gamma = 1.3

# dynamic viscosity
eta_oxygene = 42e-6
eta_propane = 117e-6 # Pa/s

zeta1 = 0
zeta2 = (1-d_pipe/d_dome)**2
zeta4 = 0
zeta5 = (1-d_pipe/d_dome)**2
zeta6 = 30
zeta7 = 2.85 #(1-d_injector/d_combustor)**2 #2.85
zeta8 = 30
zeta9 = 2.85 #2.85

density_oxygene = 1141
density_fuel = 500

bulk_modulus_oxygene = 0.2e9
bulk_modulus_fuel = 0.36e9

material_pipe = 50e9
material_dome = 50e9
material_injector = 50e9


g = 9.81 # gravity constant 
# =============================================================================
# Initial values of my variables
# =============================================================================

P0i = 40e5 # Pa
P1i = 0
P2i = 0
P3i = 40e5
P4i = 0
P5i = 0
P6i = 0
P7i = 0.01
P8i = 0
P9i = 0
m0i = 0 # kg/s
m1i = 0 # kg/s
m2i = 0
m3i = 0
m4i = 0
m5i = 0
m6i = 0
m7i = 0.0001
m8i = 0
m9i = 0.0001
m10i = 0
T10i = 300 # K
# =============================================================================
# Initialisation of my variables
# =============================================================================

P0 = P0i + density_oxygene*g*l_tank
P1 = P1i
P2 = P2i
P3 = P3i + density_fuel*g*l_tank
P4 = P4i
P5 = P5i
P6 = P6i
P7 = P7i
P8 = P8i
P9 = P9i
m0 = m0i
m1 = m1i
m2 = m2i
m3 = m3i
m4 = m4i
m5 = m5i
m6 = m6i
m7 = m7i
m8 = m8i
m9 = m9i
m10 = m10i
T10 = T10i

# =============================================================================
# Few functions for calculus
# =============================================================================

def a_sound(density, bulk_modulus,length, material):
    # if liquid, option to return something cste
    # return ( bulk_modulus / density ) / (1+(length/0.1)*(bulk_modulus/material)) # corrected speed of sound
    return  bulk_modulus / density

def Reynolds(m_flow, length, section, eta):
    m=m_flow
    if m_flow==0:
        m=1e-12
    return (m * length ) / ( section * eta )
gas1 = ct.Solution('gri30.yaml')

# Function to solve implicit problem 
def function(S):
    # /!\ to properly say it, xi is defined by the time step just before 
    # (as we solve implicit we cannot add built in function in the solver)
    # we consider the time step to be sufficiently small that :
    #    P(n+1) close to P(n) HYPOTHESIS
    lambda1 = 64/Reynolds(m1,l_pipe,A_pipe, eta_oxygene)*1/d_pipe *l_pipe
    xi1 = 0.5 * A_pipe**(-2) * ( lambda1 + zeta1 )
    
    lambda2 = 64/Reynolds(m2,l_pipe,A_pipe, eta_oxygene)*1/d_pipe *l_pipe
    xi2 = 0.5 * A_pipe**(-2) * ( lambda2 + zeta2 )
    
    lambda4 = 64/Reynolds(m4,l_pipe,A_pipe, eta_propane)*1/d_pipe *l_pipe
    xi4 = 0.5 * A_pipe**(-2) * ( lambda4 + zeta4 )
    
    lambda5 = 64/Reynolds(m5,l_pipe,A_pipe, eta_propane)*1/d_pipe *l_pipe
    xi5 = 0.5 * A_pipe**(-2) * ( lambda5 + zeta5 )
    
    lambda6 = 64/Reynolds(m6,l_dome,A_dome, eta_oxygene)*1/d_dome *l_dome
    xi6 = 0.5 * A_dome**(-2) * ( lambda6 + zeta6 )
    
    lambda8 = 64/Reynolds(m8,l_dome,A_dome, eta_propane)*1/d_dome *l_dome
    xi8 = 0.5 * A_dome**(-2) * ( lambda8 + zeta8 )
    
    lambda7 = 64/Reynolds(m7,l_injector,A_injector, eta_oxygene)*1/d_injector *l_injector
    xi7 = 0.5 * A_injector**(-2) * ( lambda7 + zeta7 )
    
    lambda9 = 64/Reynolds(m9,l_injector,A_injector, eta_propane)*1/d_injector *l_injector
    xi9 = 0.5 * A_injector**(-2) * ( lambda9 + zeta9 )
    
    gas1.Y = f'C3H8:{m9} , O2:{m7}'
    gas1.TP = T10, P7
    gas1.equilibrate('HP')
    T10n = gas1.T 
    
    return np.array([
        (S[0]-P0)/dt + S[10]*g/A_tank ,
        (S[3]-P3)/dt + S[13]*g/A_tank , 
        (S[0]-P0)/dt - a_sound(density_oxygene, bulk_modulus_oxygene,l_pipe, material_pipe) / V_pipe * (S[10] - S[11]) ,
        (S[1]-P1)/dt - a_sound(density_oxygene, bulk_modulus_oxygene,l_pipe, material_pipe) / V_pipe * (S[11] - S[12]) ,
        (S[3]-P3)/dt - a_sound(density_fuel, bulk_modulus_fuel,l_pipe, material_pipe) / V_pipe * (S[13] - S[14]) ,
        (S[4]-P4)/dt - a_sound(density_fuel, bulk_modulus_fuel,l_pipe, material_pipe) / V_pipe * (S[14] - S[15]) ,
        (S[11] - m1)/dt - A_pipe/l_pipe * (S[0] - S[1] - xi1 / density_oxygene * S[11]*np.abs(S[11])) ,
        (S[12] - m2)/dt - A_pipe/l_pipe * (S[1] - S[2] - xi2 / density_oxygene * S[12]*np.abs(S[12])) ,
        (S[14] - m4)/dt - A_pipe/l_pipe * (S[3] - S[4] - xi4 / density_fuel * S[14]*np.abs(S[14])) ,
        (S[15] - m5)/dt - A_pipe/l_pipe * (S[4] - S[5] - xi5 / density_fuel * S[15]*np.abs(S[15])) ,
        (S[2]-P2)/dt - a_sound(density_oxygene, bulk_modulus_oxygene,l_dome, material_dome) / V_dome * (S[12] - S[16]) ,
        (S[16] - m6)/dt - A_dome/l_dome * (S[2] - S[6] - xi6 / density_oxygene * S[16]*np.abs(S[16])) ,
        (S[6]-P6)/dt - a_sound(density_oxygene, bulk_modulus_oxygene,l_injector, material_injector) / V_injector * (S[16]/NI - S[17]/NI) ,
        (S[17]/NI - m7/NI)/dt - A_injector/l_injector * ( S[6] - S[7] - xi7 / density_oxygene * S[17]*np.abs(S[17])/NI**2 ) ,
        (S[5]-P5)/dt - a_sound(density_fuel, bulk_modulus_fuel,l_dome, material_dome) / V_dome * (S[15] - S[18]) ,
        (S[18] - m8)/dt - A_dome/l_dome * (S[5] - S[8] - xi8 / density_fuel * S[18]*np.abs(S[18])) ,
        (S[8]-P8)/dt - a_sound(density_fuel, bulk_modulus_fuel,l_injector, material_injector) / V_injector * (S[18]/NI - S[19]/NI) ,
        (S[19]/NI - m9/NI)/dt - A_injector/l_injector * ( S[8] - S[9] - xi9 / density_oxygene * S[19]*np.abs(S[19])/NI**2 ) ,
        S[7] - S[9] ,
        S[20] - ( np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) * A_throat )/(np.sqrt(R_cc * S[21])) * S[7] ,
        # V_combustor/( R_cc * S[21] ) * (S[7] - P7)/dt + S[20] - S[19] - S[17] - 0.5*np.sin(100*2*np.pi*t),
        # V_combustor/( R_cc * S[21] ) * (S[7] - P7)/dt + S[20] - S[19] - S[17] ,
        (S[7]- P7)/dt - R_cc * S[21] / V_combustor * ( S[17] + S[19] - S[20] + V_combustor * S[7] / R_cc * (S[21] - T10)/dt * (+1 * 1/S[21]**2) ),
        S[21] - T10n 
        # S[21] - 3000
        ])

# =============================================================================
# Plots lists
# =============================================================================
list_pressure_0 = []
list_pressure_1 = []
list_pressure_2 = []
list_pressure_3 = []
list_pressure_4 = []
list_pressure_5 = []
list_pressure_6 = []
list_pressure_7 = []
list_pressure_8 = []
list_pressure_9 = []

list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []
list_mass_flow_4 = []
list_mass_flow_5 = []
list_mass_flow_6 = []
list_mass_flow_7 = []
list_mass_flow_8 = []
list_mass_flow_9 = []
list_mass_flow_10 = []

list_temperature_10 = []

time_storage = []

start = time.time() 
print("Start computing") 

while t<0.3: 
    
    ############################################################  
    #### LOOP AND CALCULUS 
    ############################################################
    
    t = t + dt
    S_vector = [ P0, P1, P2, P3, P4, P5, P6, P7, P8, P9,
                m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10,
                T10
                ]
    root = fsolve(function, S_vector)
    P0 = root[0]
    P1 = root[1]
    P2 = root[2]
    P3 = root[3]
    P4 = root[4]
    P5 = root[5]
    P6 = root[6]
    P7 = root[7]
    P8 = root[8]
    P9 = root[9]
    m0 = root[10]
    m1 = root[11]
    m2 = root[12]
    m3 = root[13]
    m4 = root[14]
    m5 = root[15]
    m6 = root[16]
    m7 = root[17]
    m8 = root[18]
    m9 = root[19]
    m10 = root[20]
    T10 = root[21]
    
    ############################################################  
    #### SAVE 
    ############################################################
    
    list_pressure_0.append(P0)
    list_pressure_1.append(P1)
    list_pressure_2.append(P2)
    list_pressure_3.append(P3)
    list_pressure_4.append(P4)
    list_pressure_5.append(P5)
    list_pressure_6.append(P6)
    list_pressure_7.append(P7)
    list_pressure_8.append(P8)
    list_pressure_9.append(P9)
    
    list_mass_flow_0.append(m0)
    list_mass_flow_1.append(m1)
    list_mass_flow_2.append(m2)
    list_mass_flow_3.append(m3)
    list_mass_flow_4.append(m4)
    list_mass_flow_5.append(m5)
    list_mass_flow_6.append(m6)
    list_mass_flow_7.append(m7/NI)
    list_mass_flow_8.append(m8)
    list_mass_flow_9.append(m9/NI)
    list_mass_flow_10.append(m10)
    
    list_temperature_10.append(T10)
    
    time_storage.append(t)
    
end = time.time() 
print(f'time of simulation : {end - start}')   


# =============================================================================
# Plots and graphs 
# =============================================================================

plt.figure()
plt.plot(time_storage[::10], list_pressure_0[::10], marker='o', markersize=3, markevery=22)
plt.plot(time_storage[::10], list_pressure_1[::10], marker='s', markersize=3, markevery=20)
plt.plot(time_storage[::10], list_pressure_2[::10], marker='^', markersize=4, markevery=18)
plt.plot(time_storage[::10], list_pressure_6[::10], marker='*', markersize=3, markevery=15)
plt.plot(time_storage[::10], list_pressure_7[::10], marker='x', markersize=3, markevery=25)

plt.title(f'Pressure evolution in system Tank-Pipe-Pipe-Dome-Injectors-Combustor \n Oxygene Line with T changing  \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['pipe1','pipe2', 'dome', 'injector', 'combustor'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time_storage[::10], list_pressure_3[::10], marker='o', markersize=3, linestyle='-', markevery=20)
plt.plot(time_storage[::10], list_pressure_4[::10], marker='s', markersize=3, linestyle='-', markevery=22)
plt.plot(time_storage[::10], list_pressure_5[::10], marker='^', markersize=4, linestyle='-', markevery=18)
plt.plot(time_storage[::10], list_pressure_8[::10], marker='*', markersize=3, linestyle='-', markevery=15)
plt.plot(time_storage[::10], list_pressure_9[::10], marker='x', markersize=3, linestyle='-', markevery=25)


plt.title(f'Pressure evolution in system Tank-Pipe-Pipe-Dome-Injectors-Combustor \n Propane Line with T changing  \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['pipe3','pipe4', 'dome', 'injector', 'combustor'])
plt.ylabel(r'$Pressure$ (Pa)')


plt.figure()
plt.plot(time_storage[::10], list_mass_flow_0[::10], marker='o', markersize=3, markevery=28)
plt.plot(time_storage[::10], list_mass_flow_1[::10], marker='s', markersize=3, markevery=15)
plt.plot(time_storage[::10], list_mass_flow_2[::10], marker='^', markersize=3, markevery=20)
plt.plot(time_storage[::10], list_mass_flow_6[::10], marker='*', markersize=3, markevery=18)
plt.plot(time_storage[::10], list_mass_flow_7[::10], marker='x', markersize=3, markevery=20)
plt.plot(time_storage[::10], list_mass_flow_10[::10], marker='D', markersize=3, markevery=22)


plt.title(f'Mass flow evolution in system Tank-Pipe-Pipe-Dome-Injectors-Combustor \n Oxygene line with T changing   \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['tank','pipe 1','pipe 2','dome','injector','combustor out'])
plt.ylabel(r'$mass flow$ (kg/s)')
# plt.ylim(0,5)

plt.figure()
plt.plot(time_storage[::10], list_mass_flow_3[::10], marker='o', markersize=3, markevery=25)
plt.plot(time_storage[::10], list_mass_flow_4[::10], marker='s', markersize=3, markevery=18)
plt.plot(time_storage[::10], list_mass_flow_5[::10], marker='^', markersize=3, markevery=15)
plt.plot(time_storage[::10], list_mass_flow_8[::10], marker='*', markersize=3, markevery=22)
plt.plot(time_storage[::10], list_mass_flow_9[::10], marker='x', markersize=3, markevery=20)
plt.plot(time_storage[::10], list_mass_flow_10[::10], marker='D', markersize=3, markevery=28)


plt.title(f'Mass flow evolution in system Tank-Pipe-Pipe-Dome-Injectors-Combustor \n Propane Line with T changing  \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['tank','pipe 3','pipe 4','dome','injector','combustor out'])
plt.ylabel(r'$mass flow$ (kg/s)')
# plt.ylim(0,5)


print(f'Drop of pressure between injection and combustor is equal to : {int((1 - P7/P6)*100)} %')

plt.figure()
plt.plot(time_storage[::2],list_temperature_10[::2]) 
plt.title(f'Temperature evolution in the combustion chamber \n in system Tank-Pipe-Pipe-Dome-Injectors-Combustor,  with changing Temperature T10 \n')
plt.xlabel(r'Time (s)')
plt.ylabel(r'T (K)')
plt.grid()



# =============================================================================
# Fast Fourier Transformation for frequencies analysis 
# =============================================================================

from numpy.fft import fft, fftfreq
plt.figure()
# Calcul FFT
x=list_pressure_1
N= len(x)
X = fft(x)  # Transformée de fourier
freq = fftfreq(N, dt)  # Fréquences de la transformée de Fourier
# Calcul du nombre d'échantillon
N= len(x)
# On prend la valeur absolue de l'amplitude uniquement pour les fréquences positives et normalisation
X_abs = np.abs(X[:N//2])*2.0/N
# On garde uniquement les fréquences positives
freq_pos = freq[:N//2]
m=max(X_abs)
# index=X_abs.index(m)
plt.plot(freq_pos, X_abs, label="Amplitude absolue")
plt.xlim(0, 100)  # On réduit la plage des fréquences à la zone utile
# plt.ylim(0,1e5)
plt.grid()
plt.xlabel(r"Fréquence (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Transformée de Fourier")
plt.show()

# Get the indices of maximum element in numpy array
result = np.where(X_abs == np.amax(X_abs))