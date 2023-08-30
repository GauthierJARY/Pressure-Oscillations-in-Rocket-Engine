# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 15:40:34 2023

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
dt = 1e-8
time_simulation = 0.05

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
P1i = 10e5
P2i = 10e5
P3i = 40e5
P4i = 10e5
P5i = 10e5
P6i = 10e5
P7i = 10e5
P8i = 10e5
P9i = 10e5
m0i = 0.01 # kg/s
m1i = 0.01 # kg/s
m2i = 0.01
m3i = 0.01
m4i = 0.01
m5i = 0.01
m6i = 0.01
m7i = 0.01
m8i = 0.01
m9i = 0.01
m10i = 0.01
T10i = 3000 # K
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

def mu_compos(m1,m2):
    if m2 == 0 :
        mu_ox_burned, mu_ox_out, mu_fu_burned, mu_fu_out = 1,1,1,1
    else:
        MRin = m1/m2
        MRst = 2
        if MRin <= MRst:
            mu_ox_burned = MRin / (MRin + 1)
            mu_ox_out = 0
            mu_fu_burned = MRin / ( (MRin + 1) * MRst)
            mu_fu_out = (1 - MRin/MRst)/(MRin + 1)
        else :
            mu_ox_burned = MRst/ (MRin + 1)
            mu_ox_out =(MRin - MRst)/(MRin + 1)
            mu_fu_burned = 1 / (1 + MRin)
            mu_fu_out = 0
    return [mu_ox_burned, mu_ox_out, mu_fu_burned, mu_fu_out]

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
gas1 = ct.Solution('gri30.yaml')

while t<0.05: 
    
    ############################################################  
    #### CALCULUS 
    ############################################################
    
    t = t + dt
    if np.abs(t-time_simulation)<0.01:
        print('Half of simulation done \n')
        print(f'time of simulation until half of it: {time.time() - start}') 
    
    gas1.Y = f' C3H8:{ mu_compos(m7, m9)[2] * m9 } , O2:{ mu_compos(m7, m9)[0] * m7 }'
    gas1.TP = T10, P7
    gas1.equilibrate('HP')
    T10n = gas1.T
    
    P0n = P0 - dt * ( - a_sound(density_oxygene, bulk_modulus_oxygene,l_pipe, material_pipe) / V_pipe * (m0 - m1) )
    P1n = P1 - dt * ( - a_sound(density_oxygene, bulk_modulus_oxygene,l_pipe, material_pipe) / V_pipe * (m1 - m2)  ) 
    P3n = P3 - dt * ( - a_sound(density_fuel, bulk_modulus_fuel,l_pipe, material_pipe) / V_pipe * (m3 - m4) )
    P4n = P4 - dt * ( - a_sound(density_fuel, bulk_modulus_fuel,l_pipe, material_pipe) / V_pipe * (m4 - m5) )
    lambda1 = 64/Reynolds(m1,l_pipe,A_pipe, eta_oxygene)*1/d_pipe *l_pipe
    xi1 = 0.5 * A_pipe**(-2) * ( lambda1 + zeta1 )
    m1n = m1 - dt * ( - A_pipe/l_pipe * (P0 - P1 - xi1 / density_oxygene * m1*np.abs(m1) ) )
    m0n = m1n / ( 1 + ( (V_pipe / a_sound(density_oxygene, bulk_modulus_oxygene,l_pipe, material_pipe) )/ (g/A_tank) ) )
    lambda2 = 64/Reynolds(m2,l_pipe,A_pipe, eta_oxygene)*1/d_pipe *l_pipe
    xi2 = 0.5 * A_pipe**(-2) * ( lambda2 + zeta2 )
    m2n = m2 - dt * ( - A_pipe/l_pipe * (P1 - P2 - xi2 / density_oxygene * m2*np.abs(m2) ) )
    lambda4 = 64/Reynolds(m4,l_pipe,A_pipe, eta_propane)*1/d_pipe *l_pipe
    xi4 = 0.5 * A_pipe**(-2) * ( lambda4 + zeta4 )
    m4n = m4 - dt * ( - A_pipe/l_pipe * (P3 - P4 - xi4 / density_fuel * m4*np.abs(m4)) )
    m3n = m4n / ( 1 + ( (V_pipe / a_sound(density_fuel, bulk_modulus_fuel,l_pipe, material_pipe) )/ (g/A_tank) ) )
    lambda5 = 64/Reynolds(m5,l_pipe,A_pipe, eta_propane)*1/d_pipe *l_pipe
    xi5 = 0.5 * A_pipe**(-2) * ( lambda5 + zeta5 )
    m5n = m5 - dt * ( - A_pipe/l_pipe * (P4 - P5 - xi5 / density_fuel * m5*np.abs(m5)) )
    P2n = P2 - dt * ( - a_sound(density_oxygene, bulk_modulus_oxygene,l_dome, material_dome) / V_dome * (m2 - m6) )
    lambda6 = 64/Reynolds(m6,l_dome,A_dome, eta_oxygene)*1/d_dome *l_dome
    xi6 = 0.5 * A_dome**(-2) * ( lambda6 + zeta6 )
    m6n = m6 - dt * ( - A_dome/l_dome * (P2 - P6 - xi6 / density_oxygene * m6*np.abs(m6)) )
    P6n = P6 - dt * ( - a_sound(density_oxygene, bulk_modulus_oxygene,l_injector, material_injector) / V_injector * (m6/NI - m7/NI) )
    lambda7 = 64/Reynolds(m7,l_injector,A_injector, eta_oxygene)*1/d_injector *l_injector
    xi7 = 0.5 * A_injector**(-2) * ( lambda7 + zeta7 )
    m7n = m7 - dt * NI * ( - A_injector/l_injector * ( P6 - P7 - xi7 / density_oxygene * m7*np.abs(m7)/NI**2 ) )
    P5n = P5 - dt * ( - a_sound(density_fuel, bulk_modulus_fuel,l_dome, material_dome) / V_dome * (m5 - m8) )
    lambda8 = 64/Reynolds(m8,l_dome,A_dome, eta_propane)*1/d_dome *l_dome
    xi8 = 0.5 * A_dome**(-2) * ( lambda8 + zeta8 )
    m8n = m8 - dt * (- A_dome/l_dome * (P5 - P8 - xi8 / density_fuel * m8*np.abs(m8)) )
    P8n = P8 - dt * ( - a_sound(density_fuel, bulk_modulus_fuel,l_injector, material_injector) / V_injector * (m8/NI - m9/NI) )
    lambda9 = 64/Reynolds(m9,l_injector,A_injector, eta_propane)*1/d_injector *l_injector
    xi9 = 0.5 * A_injector**(-2) * ( lambda9 + zeta9 )
    m9n = m9 - dt * NI * ( - A_injector/l_injector * ( P8 - P9 - xi9 / density_oxygene * m9*np.abs(m9)/NI**2 ) ) 
    T10n = 3000
    # P7n = P7 - dt * (m10 - m9 - m7) * R_cc * T10 / V_combustor
    P7n = P7 + dt * R_cc * T10 / V_combustor * ( m7 + m9 - m10 + V_combustor * P7 / R_cc * (T10n - T10)/dt * (+1 * 1/T10**2) )
    P9n = P7n
    m10n = - ( - ( np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) * A_throat )/(np.sqrt(R_cc * T10n)) * P7n ) 

    ############################################################  
    #### LOOP
    ############################################################  
    P0 = P0n
    P1 = P1n
    P2 = P2n
    P3 = P3n
    P4 = P4n
    P5 = P5n
    P6 = P6n
    P7 = P7n
    P8 = P8n
    P9 = P9n
    m0 = m0n
    m1 = m1n
    m2 = m2n
    m3 = m3n
    m4 = m4n
    m5 = m5n
    m6 = m6n
    m7 = m7n
    m8 = m8n
    m9 = m9n
    m10 = m10n
    T10 = T10n
    
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
plt.plot(time_storage[::10],list_pressure_0[::10])
plt.plot(time_storage[::10],list_pressure_1[::10])
plt.plot(time_storage[::10],list_pressure_2[::10])
plt.plot(time_storage[::10],list_pressure_6[::10])
plt.plot(time_storage[::10],list_pressure_7[::10])

plt.title(f'Pressure evolution Tank-Pipe-Pipe-Dome-Injectors-Combustor \n Oxygene Line with T3 = {T3}K \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['pipe1','pipe2', 'dome', 'injector', 'combustor'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time_storage[::10],list_pressure_3[::10])
plt.plot(time_storage[::10],list_pressure_4[::10])
plt.plot(time_storage[::10],list_pressure_5[::10])
plt.plot(time_storage[::10],list_pressure_8[::10])
plt.plot(time_storage[::10],list_pressure_9[::10])

plt.title(f'Pressure evolution for Tank-Pipe-Pipe-Dome-Injectors-Combustor \n Propane Line with T3 = {T3}K \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['pipe1','pipe2', 'dome', 'injector', 'combustor'])
plt.ylabel(r'$Pressure$ (Pa)')


plt.figure()
plt.plot(time_storage[::10],list_mass_flow_0[::10])
plt.plot(time_storage[::10],list_mass_flow_1[::10])
plt.plot(time_storage[::10],list_mass_flow_2[::10])
plt.plot(time_storage[::10],list_mass_flow_6[::10])
plt.plot(time_storage[::10],list_mass_flow_7[::10])
plt.plot(time_storage[::10],list_mass_flow_10[::10])

plt.title(f'massflow evolution Tank-Pipe-Pipe-Dome-Injectors-Combustor \n Oxygene with T3 = {T3}K  \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['tank','pipe','pipe','dome','injector','combustor out'])
plt.ylabel(r'$mass flow$ (kg/s)')
plt.ylim(0,5)

plt.figure()
plt.plot(time_storage[::10],list_mass_flow_3[::10])
plt.plot(time_storage[::10],list_mass_flow_4[::10])
plt.plot(time_storage[::10],list_mass_flow_5[::10])
plt.plot(time_storage[::10],list_mass_flow_8[::10])
plt.plot(time_storage[::10],list_mass_flow_9[::10])
plt.plot(time_storage[::10],list_mass_flow_10[::10])

plt.title(f'massflow evolution Tank-Pipe-Pipe-Dome-Injectors-Combustor \n Propane Line with T3 = {T3}K  \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['tank','pipe','pipe','dome','injector','combustor out'])
plt.ylabel(r'$mass flow$ (kg/s)')
plt.ylim(0,5)


print(f'Drop of pressure between injection and combustor is equal to : {int((1 - P7/P6)*100)} %')

plt.figure()
plt.plot(time_storage[::2],list_temperature_10[::2]) 
plt.title(f'temperature evolution in Combustor fed at changing temperature \n')
plt.xlabel(r'$time$ (s)')
plt.ylabel(r'$T$ (K)')

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
plt.xlim(0, 500)  # On réduit la plage des fréquences à la zone utile
plt.ylim(0,1e5)
plt.grid()
plt.xlabel(r"Fréquence (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Transformée de Fourier")
plt.show()

# Get the indices of maximum element in numpy array
result = np.where(X_abs == np.amax(X_abs))