# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:32:54 2023

@author: gauth
"""

# =============================================================================
# READ ME 
#    Important code : 1st time Lumped Parameter Method was found in this internship and implemented
#    it is the 1st try once i discovered the lumped hyptohesis 
#    very close to A1 TPPT Chiara  
# =============================================================================

# =============================================================================
# Imported Module
# =============================================================================

import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys

# =============================================================================
# Geometric Input
# =============================================================================

l_pipe_1 = 2 #meter
l_pipe_2 = 2
l_tank_1 = 10
l_tank_2 = 10

d_tank_1 = 1 #meter 
d_tank_2 = 1
d_pipe_1 = 0.5
d_pipe_2 = 0.5

A_pipe_1 = (d_pipe_1/2)**2 * 3.14 #m^2
A_pipe_2 = (d_pipe_2/2)**2 * 3.14
A_tank_1 = (d_tank_1/2)**2 * 3.14
A_tank_2 = (d_tank_2/2)**2 * 3.14

V_tank_1 = l_tank_1 * A_tank_1 #m^3
V_tank_2 = l_tank_2 * A_tank_2
V_pipe_1 = l_pipe_1 * A_pipe_1
V_pipe_2 = l_pipe_2 * A_pipe_2


# =============================================================================
# Physical constants Input
# =============================================================================

type_simulation = 1 # 0 : gas & 1 : liquid

rho0i = 1197 #kg/m^3, water density 
rho1i = 1197
rho2i = 1197

t = 0 #s
dt = 0.000001 # time step of the simulation explicite Euler

gamma = 1 # isentropic coefficient
nu = 1e-5 # cinematic viscosity water
eta = 1e-3 # dynamic viscosity water

zeta_geo_1 = 0.5625 # loss coefficient from report 
zeta_geo_2 = 9 # loss coefficient from report

material_pipe = 50e9 #Pa , Young Modulus aluminium
beta = 2.2e9 #Pa, bulk modulus water

# =============================================================================
# Initialisation of the variable
# =============================================================================

# initial mass flow and pressure

P0i = 0.6e6 #Pa
P1i = 0.4e6
P2i = 0.5e6

m0i = 0 #kg/s
m1i = 0
m2i = 0
m3i = 0

# variables 

m0_prev = m0i
m1_prev = m1i
m2_prev= m2i
m3_prev = m3i

P0 = P0i + rho0i*10*l_tank_1
P1 = P1i
P2 = P2i + rho2i*10*l_tank_2

# =============================================================================
# List for saving evolution
# =============================================================================

list_pressure_0 = []
list_pressure_1 = []
list_pressure_2 = []

list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []

time = []

# =============================================================================
# Functions 
# =============================================================================

def a_sound(Pressure_current,Pressure_init,rho_init):
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        return ( beta / rho_init ) / (1+(d_pipe_1/0.02)*(beta/material_pipe))
        
    
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
# Loop for computing 
# =============================================================================


while t<0.05:
    
    t = t + dt

    # =========================================================================
    #     Loop to calculus
    # =========================================================================
    
    m0_prev = m1_prev /( 1 + ( (V_pipe_1/a_sound(P1, P1i, rho1i))/(9.81/A_tank_1) ) )
    m0_next = m0_prev
    # P0n = P0 + dt * ( 9.81 / A_tank_1 ) * (- m0_prev)
    P0n = P0 + dt * ( (m0_prev-m1_prev) / (V_pipe_1 ) * a_sound(P1,P1i,rho1i) )
    P1n = P1 + dt * ( (m1_prev-m2_prev) / (V_pipe_1 ) * a_sound(P1,P1i,rho1i) )
    P2n = P2 + dt * ( ( (m2_prev) * 9.81 / A_tank_2)  )

    zeta_fluid_1 = 64/Reynolds(m1_prev,l_pipe_1,A_pipe_1)*1/d_pipe_1 *l_pipe_1
    xi1 = 0.5 * A_pipe_1**(-2) * ( zeta_fluid_1 + zeta_geo_1 )
    m1_next = m1_prev + dt * ( (A_pipe_1/l_pipe_1) * ( P0 - P1 - ( xi1 / rho(P1,P1i,rho1i) ) * m1_prev * np.abs(m1_prev) ) ) 
    
    zeta_fluid_2 = 64/Reynolds(m2_prev,l_pipe_2,A_pipe_2)*2/d_pipe_2 *l_pipe_2
    xi2 = 0.5 * A_pipe_2**(-2) * ( zeta_fluid_2 + zeta_geo_2 )
    m2_next = m2_prev + dt * ( (A_pipe_1/l_pipe_2) * ( P1 - P2 - ( xi2 / rho(P1,P1i,rho1i) ) * m2_prev * np.abs(m2_prev) ) ) 
    
    # =========================================================================
    #     Loop to integrate
    # =========================================================================
    
    P0 = P0n
    P1 = P1n
    P2 = P2n

    m0_prev = m0_next
    m1_prev = m1_next
    m2_prev = m2_next

    # =========================================================================
    #     Loop to save
    # =========================================================================
    
    list_pressure_0.append(P0n)
    list_pressure_1.append(P1)
    list_pressure_2.append(P2n)

    list_mass_flow_0.append(m0_prev)
    list_mass_flow_1.append(m1_prev)
    list_mass_flow_2.append(m2_prev)
    
    time.append(t)

    # =========================================================================
    #    Eigen Values Analysis 
    # =========================================================================
    
    if np.abs(t-0.1)<0.5*dt:
        
        celerity_sound = np.sqrt(a_sound(P1,P1i,rho1i))
        print(f'time : {t} \n P0 : {P0}\n P1 : {P1} \n P2 : {P2}\n m0 : {m0_prev} \n m1 : {m1_prev}\n m2 : {m2_prev} \n a : {celerity_sound}')

        celerity_sound = np.sqrt(a_sound(P1,P1i,rho1i))
        print(f'a speed sound non corrected : {celerity_sound*np.sqrt(1+(l_pipe_1/0.1)*(2.2e9/50e9))}')
        print(f'a speed sound corrected : {celerity_sound}')
        V1=V_pipe_1 # we proved it in our report 
        lin1 = A_pipe_1/l_pipe_1*xi1/rho1i*np.abs(m1_prev)
        lin2 = A_pipe_2/l_pipe_2*xi2/rho2i*np.abs(m2_prev)        
        # lin1 = A_pipe_1/l_pipe_1*xi1/rho1i * 1/1e-6
        # lin2 = A_pipe_2/l_pipe_2*xi2/rho2i * 1/1e-6

        A = np.array( [ [0,0,0,-9.81/A_tank_1,0],[0,0,0,celerity_sound**2/V1,-celerity_sound**2/V1],[0,0,0,0,9.81/A_tank_2],[A_pipe_1/l_pipe_1,-A_pipe_1/l_pipe_1,0,lin1,0],[0,A_pipe_2/l_pipe_2,-A_pipe_2/l_pipe_2,0,lin2]])
        eigen_values_matrix = np.linalg.eigvals(-A)
        print(f'eigen values are : {eigen_values_matrix}')
        eigen_frequencies = np.imag(eigen_values_matrix[0])/(2*3.14)
        print(f'eigen frequencies are : {eigen_frequencies} Hz')
        print(f'harmonic frequencies : {celerity_sound/(2*2*l_pipe_1)}') # theory of wind tubes

#end of the loop


# =============================================================================
# Plots
# =============================================================================

plt.figure()
plt.plot(time[::10], list_pressure_0[::10], marker='^', markersize=4, markevery=800)
plt.plot(time[::10], list_pressure_1[::10], marker='x', markersize=4, markevery=200)
plt.plot(time[::10], list_pressure_2[::10], marker='*', markersize=4, markevery=800)

plt.title(f'Pressure Evolution in the Tanks and the Pipelines, System TPPT, filled with water\n')
plt.xlabel(r'Time (s)')
plt.ylabel(r'Pressure (Pa)')

# Placing the legend in the bottom right corner
plt.legend(['Tank n°1 Pressure P0', 'Pipe n°2 Pressure P1', 'Tank n°2 Pressure P2'], loc='lower right')

plt.show()

plt.figure()
plt.plot(time[::10],list_mass_flow_0[::10], marker='^', markersize=4,markevery=160)
plt.plot(time[::10],list_mass_flow_1[::10], marker='x', markersize=4,markevery=200)
plt.plot(time[::10],list_mass_flow_2[::10], marker='*', markersize=4, markevery=200)

plt.title(f'Mass flow Evolution in the Tanks and the Pipelines, System TPPT, filled with water\n')
plt.xlabel(r'Time (s)')
plt.ylabel(r'Mass flow (kg/s)')
plt.legend(['Massflow in Tank n°1 $\dot{m_0}$','Massflow in Pipe n°1 $\dot{m_1}$','Massflow in Pipe n°2 $\dot{m_2}$'])


# =============================================================================
# Frequency Analysis : FFT 
# =============================================================================

from numpy.fft import fft, fftfreq

plt.figure()

# Calcul FFT

x=list_pressure_1
x=x+ [0]*int(len(x)/8) # zero padding 
N = len(x)

X_amplitude = fft(x)  # Transformée de fourier
frequencies_of_FFT = fftfreq(N, dt)  # Fréquences de la transformée de Fourier
# Calcul du nombre d'échantillon
# On prend la valeur absolue de l'amplitude uniquement pour les fréquences positives et normalisation
X__amplitude_absolute_normalized = np.abs(X_amplitude[:N//2])*2.0/N
# On garde uniquement les fréquences positives
frequencies_of_FFT_positive = frequencies_of_FFT[:N//2]


plt.plot(frequencies_of_FFT_positive, X__amplitude_absolute_normalized, label="Amplitude absolue")
plt.ylim(0,0.2e6)
plt.xlim(0, 200)  # On réduit la plage des fréquences à la zone utile
plt.grid()
plt.xlabel(r"Frequency (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Fast Fourier Transformation of the pressure $P_1$")
plt.show()



