# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 15:07:48 2023

@author: gauth
"""

# =============================================================================
# README
#      comparable with chiara results in example A2
#     check for equations and boundary conditions
# =============================================================================


# =============================================================================
# Imported Modules
# =============================================================================
import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys

# =============================================================================
# Geometric Input 
# =============================================================================
l_pipe_1 = 4
l_tank_1 = 1
l_tank_2 = 1

d_tank_1 = 1
d_tank_2 = 1
d_pipe_1 = 0.5

A_pipe_1 = (d_pipe_1/2)**2 * 3.14
A_tank_1 = (d_tank_1/2)**2 * 3.14
A_tank_2 = (d_tank_2/2)**2 * 3.14

V_tank_1 = l_tank_1 * A_tank_1
V_tank_2 = l_tank_2 * A_tank_2
V_pipe_1 = l_pipe_1 * A_pipe_1

# =============================================================================
# Physical Constants Input
# =============================================================================
rho0i = 1197
rho1i = 1197
rho2i = 1197

t = 0
dt = 0.0001
type_simulation = 1 # 0 : gas & 1 : liquid

gamma = 1 # isentropic coefficient
nu = 1e-5 # cinematic viscosity water
eta = 1e-3 # dynamic viscosity water

zeta_geo_1 = 0.5625 # loss coeff due to geometry
zeta_geo_2 = 9

material_pipe = 50e9 #Pa , Young Modulus aluminium
beta = 2.2e9 #Pa, bulk modulus water

# =============================================================================
# Initialisation 
# =============================================================================
P0i = 0.2e6 #Pa 
P1i = 0.1e6
P2i = 0.1e6

m1i = 0 #kg/s
m1_prev = m1i

P0 = P0i + rho0i*9.81*l_tank_1
P1 = P1i
P2 = P2i 

# =============================================================================
# Lists to save 
# =============================================================================
list_pressure_0 = []
list_pressure_1 = []
list_pressure_2 = []
list_massflow_0 = []
time = []

# =============================================================================
# Functions
# =============================================================================

def a_sound(Pressure_current,Pressure_init,rho_init,diameter):
    # if liquid, option to return something cste
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        # corrected value of sound with the confinement effect
        return ( beta / rho_init ) / (1+(diameter/0.02)*(beta/material_pipe))      

def rho(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 : 
        return rho_init * (Pressure_current/Pressure_init)**(1/gamma)
    else :
        return rho_init

def Reynolds(m_flow, length, section):
    m=m_flow
    if m_flow==0:
        m=1e-12
    return (m * length ) / ( section * eta )

# =============================================================================
# LOOP
# =============================================================================
while t<20:
    t = t + dt  
    # =========================================================================
    #     # loop to calculus
    # =========================================================================
    P0n = P0 + dt * ( 9.81 / A_tank_1 ) * (- m1_prev)
    P2n = P2 + dt * ( ( (m1_prev) * 9.81 / A_tank_2)  )
    P1n = P1 + dt * ( (m1_prev ) / (V_pipe_1) * a_sound(P1,P1i,rho1i,d_pipe_1) )

    zeta_fluid_1 = 64/Reynolds(m1_prev,l_pipe_1,A_pipe_1)*1/d_pipe_1 *l_pipe_1
    xi1 = 0.5 * A_pipe_1**(-2) * ( zeta_fluid_1 + zeta_geo_1 )
    m1_next = m1_prev + dt * ( (A_pipe_1/l_pipe_1) * ( P0 - P2 - ( xi1 / rho(P2,P2i,rho2i) ) * m1_prev * np.abs(m1_prev) ) ) 
    sound_celerity = np.sqrt(a_sound(P2,P2i,rho2i,d_pipe_1))

    # =========================================================================
    #     # loop to integrate 
    # =========================================================================
    P0 = P0n
    P2 = P2n
    m1_prev = m1_next
    
    # =========================================================================
    #     # loop to save 
    # =========================================================================
    list_pressure_0.append(P0n)
    list_pressure_1.append(P1n)
    list_pressure_2.append(P2n)
    list_massflow_0.append(m1_prev)
    time.append(t)
    
# =============================================================================
#     Plots
# =============================================================================

plt.figure()
plt.plot(time[::1000],list_pressure_0[::1000],marker='o', markevery=10,markersize=4)
plt.plot(time[::1000],list_pressure_2[::1000],marker='^', markevery=10,markersize=4)
plt.title(f'Pressure evolution Tank-Pipe-Tank system, Filled with water\n')
plt.xlabel(r'Time (s)')
plt.legend(['Tank n°1 pressure $P_0$', ' Tank n°2 pressure $P_2$'])
plt.ylabel(r'Pressure (Pa)')

plt.figure()
plt.plot(time[::1000],list_massflow_0[::1000], linestyle='dashed')

plt.title(f'Mass flow evolution Tank-Pipe-Tank system, Filled with water\n')
plt.xlabel(r'Time (s)')
plt.legend(['Mass flow in the pipe'])
plt.ylabel(r'Mass flow (kg/s)')
plt.minorticks_on()


# =============================================================================
# FFT Analysis 
# =============================================================================
from numpy.fft import fft, fftfreq
import matplotlib.pyplot as plt

# Liste de listes
pressure_lists = [list_pressure_0, list_massflow_0, list_pressure_2]  # Ajoutez vos autres listes ici

plt.figure()
for i, x in enumerate(pressure_lists):
    # Calcul FFT
    x_extended = x + [0] * int(0.02 * len(x))
    N = len(x_extended)
    X = fft(x_extended)  # Transformée de fourier
    freq = fftfreq(N, dt)  # Fréquences de la transformée de Fourier
    # Calcul du nombre d'échantillon
    N = len(x_extended)
    # On prend la valeur absolue de l'amplitude uniquement pour les fréquences positives et normalisation
    X_abs = np.abs(X[:N//2]) * 2.0 / N
    # On garde uniquement les fréquences positives
    freq_pos = freq[:N//2]
    # Affichage des graphiques FFT sur le même graphique
    plt.semilogy(freq_pos, X_abs)

plt.xlim(0.001, 2)  # On réduit la plage des fréquences à la zone utile
plt.ylim(10, 50000)  # On réduit la plage des fréquences à la zone utile
plt.grid()
plt.xlabel(r"Frequency (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("FFT of the signals in modelling Tank-Pipe-Tank")
plt.legend( ["Pressure tank n°1","Mass flow pipe","Pressure tank n°2"])
plt.tight_layout()
plt.show()

