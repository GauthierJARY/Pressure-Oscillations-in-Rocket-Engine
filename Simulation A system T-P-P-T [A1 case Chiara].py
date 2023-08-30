# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 10:57:02 2023

@author: gauth
"""

# =============================================================================
# READ ME 
#     2 pipes implemented with chiara equations 
#     eigen values analysis
#     exactly the model of chiara in her book
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
dt = 0.0000005 # time step of the simulation explicite Euler

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
# Loop of computing
# =============================================================================

while t<0.2:
    
    t = t + dt
     
    P0n = P0 + dt * ( 9.81 / A_tank_1 ) * (m0_prev - m1_prev)
    P1n = P1 + dt * ( (m1_prev-m2_prev ) / (V_pipe_1) * a_sound(P1,P1i,rho1i) )
    P2n = P2 + dt * ( ( (m2_prev- m3_prev) * 9.81 / A_tank_2)  )

    zeta_fluid_1 = 64/Reynolds(m1_prev,l_pipe_1,A_pipe_1)*1/d_pipe_1 *l_pipe_1
    xi1 = 0.5 * A_pipe_1**(-2) * ( zeta_fluid_1 + zeta_geo_1 )
    m1_next = m1_prev + dt * ( (A_pipe_1/l_pipe_1) * ( P0 - P1 - ( xi1 / rho(P1,P1i,rho1i) ) * m1_prev * np.abs(m1_prev) ) ) 
    
    zeta_fluid_2 = 64/Reynolds(m2_prev,l_pipe_2,A_pipe_2)*2/d_pipe_2 *l_pipe_2
    xi2 = 0.5 * A_pipe_2**(-2) * ( zeta_fluid_2 + zeta_geo_2 )
    m2_next = m2_prev+ dt * ( (A_pipe_1/l_pipe_2) * ( P1 - P2 - ( xi2 / rho(P1,P1i,rho1i) ) * m2_prev* np.abs(m2_prev) ) ) 

    m0_next = m0_prev + dt * ( 1/l_tank_1 * (-A_pipe_1*P0) - rho(P0,P0i,rho0i) * 10 * A_tank_1)    
    
    # loop to integrate 
    P0 = P0n
    P1 = P1n
    P2 = P2n
    m1_prev = m1_next
    m2_prev= m2_next
    # m0_prev = m0_next
    # m3_prev = m3_next
    
    # loop to save 
    list_pressure_0.append(P0n)
    list_pressure_1.append(P1n)
    list_pressure_2.append(P2n)
    list_mass_flow_0.append(m0_prev)
    list_mass_flow_1.append(m1_prev)
    list_mass_flow_2.append(m2_prev)
    list_mass_flow_3.append(m3_prev)
    time.append(t)
    
    # Analysis 
    
    if np.abs(t-0.1)<0.5*dt:
        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        
        print(f'time : {t} \n P0 : {P0}\n P1 : {P1} \n P2 : {P2}\n m0 : {m0_prev} \n m1 : {m1_prev}\n m2 : {m2_prev} \n a : {aa}')

        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        print(f'a speed sound non corrected : {aa*np.sqrt(1+(l_pipe_1/0.1)*(2.2e9/50e9))}')
        print(f'a speed sound corrected : {aa}')
        V1=V_pipe_1
        lin1 = A_pipe_1/l_pipe_1*xi1/rho1i*np.abs(m1_prev)
        # lin1 = A_pipe_1/l_pipe_1*xi1/rho1i * 1/1e-6
        lin2 = A_pipe_2/l_pipe_2*xi2/rho2i*np.abs(m2_prev)
        # lin2 = A_pipe_2/l_pipe_2*xi2/rho2i * 1/1e-6
        A = np.array( [ [0,0,0,-9.81/A_tank_1,0],[0,0,0,aa**2/V1,-aa**2/V1],[0,0,0,0,9.81/A_tank_2],[A_pipe_1/l_pipe_1,-A_pipe_1/l_pipe_1,0,lin1,0],[0,A_pipe_2/l_pipe_2,-A_pipe_2/l_pipe_2,0,lin2]])
        vap = np.linalg.eigvals(A)
        print(f'eigen values are : {vap}')
        freqq = np.imag(vap[0])/(2*3.14)
        print(f'eigen frequencies are : {freqq} Hz')
        print(f'harmonic frequencies : {aa/(2*2*l_pipe_1)}')


# =============================================================================
# Plots
# =============================================================================

plt.figure()
plt.plot(time[::1000],list_pressure_0[::1000])
plt.plot(time[::1000],list_pressure_1[::1000])
plt.plot(time[::1000],list_pressure_2[::1000])
plt.title(f'Liquid Oxygen model : Tank Pipe Tank \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['Tank n°1 pressure P0','Pipe pressure P1', ' Tank n°2 pressure P2'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time[::1000],list_mass_flow_0[::1000])
plt.plot(time[::1000],list_mass_flow_1[::1000])
plt.plot(time[::1000],list_mass_flow_2[::1000])
plt.plot(time[::1000],list_mass_flow_3[::1000])

plt.title(f'mass evolution in filling of Liquid Oxygen of tanks \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m0','m1_prev','m2','m3'])
plt.ylabel(r'$mass flow$ (kg/s)')


# =============================================================================
# Frequency Analysis : FFT Computing 
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
plt.xlim(0, 200)  # On réduit la plage des fréquences à la zone utile
plt.ylim(0,500000)
plt.grid()
plt.xlabel(r"Frequency (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Fast Fourier Transformation")
plt.show()

