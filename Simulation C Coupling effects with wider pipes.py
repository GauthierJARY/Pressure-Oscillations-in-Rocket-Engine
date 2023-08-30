# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:39:43 2023

@author: gauth
"""


# =============================================================================
# README
#     # interesting code for report with the try to implement a cavity and see eigen values related to it 
#     be careful with singular loss charge coefficient 
#     dome is like a big pipe, interesting ! done later in the internship
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

l_pipe_1 = 2
d_pipe_1 = 0.5
A_pipe_1 = (d_pipe_1/2)**2 * 3.14
V_pipe_1 = l_pipe_1 * A_pipe_1

l_pipe_2 = 2
d_pipe_2 = 0.5
A_pipe_2 = (d_pipe_2/2)**2 * 3.14
V_pipe_2 = l_pipe_2 * A_pipe_2

l_cavity_1 = 2
d_cavity_1 = 1
A_cavity_1 = (d_cavity_1/2)**2 * 3.14
V_cavity_1 = l_cavity_1 * A_cavity_1

l_cavity_2 = 2
d_cavity_2 = 1
A_cavity_2 = (d_cavity_2/2)**2 * 3.14
V_cavity_2 = l_cavity_2 * A_cavity_2

l_tank_1 = 10
d_tank_1 = 1
A_tank_1 = (d_tank_1/2)**2 * 3.14
V_tank_1 = l_tank_1 * A_tank_1


# =============================================================================
# Physical constants Input
# =============================================================================

type_simulation = 1 # 0 : gas & 1 : liquid

rho0i = 1197
rho1i = 1197
rho2i = 1197
rho3i = 1197

t = 0
dt = 0.0000005

gamma = 1 # isentropic coefficient
nu = 1e-5 # cinematic viscosity water
eta = 1e-3 # dynamic viscosity water

zeta_geo_0 = 1 # loss coefficient
zeta_geo_1 = 0.0
zeta_geo_2 = 0.5625
zeta_geo_3 = 0

material_pipe = 50e9 #Pa , Young Modulus aluminium
beta = 2.2e9 #Pa, bulk modulus water

# =============================================================================
# Lists to save 
# =============================================================================

list_pressure_0 = [] 
list_pressure_1 = []
list_pressure_2 = [] 
list_pressure_3 = [] 
list_pressure_ext = []

list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []

time = []

# =============================================================================
# Initialisation 
# =============================================================================

P0i = 0.6e6
P1i = 0.4e6
P2i = 0.4e6
P3i = 0.4e6
Patmi = 0.1e6

m0i = 0
m1i = 0
m2i = 0
m3i = 0

m0_prev = m0i
m1_prev = m1i
m2_prev = m2i
m3_prev = m3i

P0 = P0i + 9.81 * rho0i* l_tank_1
P1 = P1i
P2 = P2i
P3 = P3i

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

while t<0.05:
    
    
    t = t + dt

    P0n = P0 + dt * ( 9.81 / A_tank_1 ) * (- m0_prev)
    P1n = P1 + dt * ( (m0_prev-m1_prev) / (V_pipe_1) * a_sound(P1,P1i,rho1i,d_pipe_1) )
    P2n = P2 + dt * ( (m1_prev-m2_prev) / (V_pipe_2) * a_sound(P2,P2i,rho2i,d_pipe_1) )
    P3n = P3 + dt * ( (m2_prev-m3_prev) / (V_cavity_2) * a_sound(P2,P2i,rho1i,d_cavity_1) )
    Patmn = Patmi
    
    zeta_fluid_0 = 64/Reynolds(m0_prev,l_pipe_1,A_pipe_1)*1/d_pipe_1 *l_pipe_1
    xi0 = 0.5 * A_pipe_1**(-2) * ( zeta_fluid_0 + zeta_geo_0 )
    m0_next = m0_prev + dt * ( (A_pipe_1/l_pipe_1) * ( P0 - P1 - ( xi0 / rho(P0,P0i,rho0i) ) * m0_prev * np.abs(m0_prev) ) ) 
    
    zeta_fluid_1 = 64/Reynolds(m1_prev,l_pipe_2,A_pipe_2)*1/d_pipe_2 *l_pipe_2
    xi1 = 0.5 * A_pipe_2**(-2) * ( zeta_fluid_1 + zeta_geo_1 )
    m1_next = m1_prev + dt * ( (A_pipe_2/l_pipe_2) * ( P1 - P2 - ( xi1 / rho(P1,P1i,rho1i) ) * m1_prev * np.abs(m1_prev) ) ) 
    
    zeta_fluid_2 = 64/Reynolds(m2_prev,l_cavity_2,A_cavity_2)*1/d_cavity_2 *l_cavity_2
    xi2 = 0.5 * A_pipe_2**(-2) * ( zeta_fluid_2 + zeta_geo_2 )
    m2_next = m2_prev + dt * ( (A_cavity_2/l_cavity_2) * ( -P3 + P2 - ( xi2 / rho(P3,P3i,rho3i) ) * m2_prev * np.abs(m2_prev) ) ) 
    
    zeta_fluid_3 = 64/Reynolds(m2_prev,l_cavity_2,A_cavity_2)*1/d_cavity_2 *l_cavity_2
    xi3 = 0.5 * A_pipe_2**(-2) * ( zeta_fluid_3 + zeta_geo_2 )
    m3_next = m3_prev + dt * ( (A_cavity_2/l_cavity_2) * ( -Patmi + P3 - ( xi3 / rho(Patmi,Patmi,rho0i) ) * m3_prev * np.abs(m3_prev) ) ) 
    

    # =========================================================================
    # Integrate
    # =========================================================================

    P0 = P0n
    P1 = P1n
    P2 = P2n
    P3 = P3n
    m0_prev = m0_next
    m1_prev = m1_next
    m2_prev = m2_next
    m3_prev = m3_next
    Patmi = Patmn
    # =========================================================================
    # Save
    # =========================================================================
    list_pressure_0.append(P0) 
    list_pressure_1.append(P1)
    list_pressure_2.append(P2)  
    list_pressure_3.append(P3) 
    list_pressure_ext.append(Patmn)
    list_mass_flow_0.append(m0_prev)
    list_mass_flow_1.append(m1_prev)
    list_mass_flow_2.append(m2_prev)
    list_mass_flow_3.append(m3_prev)
    time.append(t)
    
    # =========================================================================
    #     Analyse
    # =========================================================================
    if np.abs(t-0.05)<0.5*dt:
        aa = np.sqrt(a_sound(P1,P1i,rho1i,d_pipe_1))
        print(f'time : {t} \n P0 : {P0} \n P1 : {P1} \n P2 : {P2}\n P3 : {P3} \n P_extérieur : {Patmi} \n m0 : {m0_prev}\n m1 : {m1_prev}\n m2 : {m2_prev} \n m3 : {m3_prev}\n a : {aa}')

        aa = np.sqrt(a_sound(P3,P3i,rho3i,d_pipe_1))
        print(f'a speed sound non corrected : {aa*np.sqrt(1+(d_pipe_1/0.02)*(2.2e9/50e9))}')
        print(f'a speed sound corrected : {aa}')
        aa = np.sqrt(a_sound(P3,P3i,rho3i,d_cavity_1))
        print(f'a speed sound non corrected : {aa*np.sqrt(1+(d_cavity_1/0.02)*(2.2e9/50e9))}')
        print(f'a speed sound corrected cavity: {aa}')
# end of loop computing 



# =============================================================================
# Plots
# =============================================================================

plt.figure()
plt.plot(time[::1000],list_pressure_0[::1000])
plt.plot(time[::1000],list_pressure_1[::1000])
plt.plot(time[::1000],list_pressure_2[::1000])
plt.plot(time[::1000],list_pressure_3[::1000])
plt.plot(time[::1000],list_pressure_ext[::1000])
plt.title(f'Pressure evolution in system T-P-P-D-D, Filled with water, with outlet on atmosphere \n')
plt.xlabel(r'Time (s)')
plt.legend(['Tank $P_0$','Pipe n°1 $P_1$','Pipe n°2 $P_2$','Pipe n°3 $P_3$', 'Pipe n°4 $P_4 = P_{atm}$'], fontsize=7,loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
plt.ylabel(r'Pressure (Pa)')

plt.figure()
plt.plot(time[::100],list_mass_flow_0[::100])
plt.plot(time[::100],list_mass_flow_1[::100])
plt.plot(time[::100],list_mass_flow_2[::100])
plt.plot(time[::100],list_mass_flow_3[::100])
plt.title(f'Mass flow evolution in system T-P-P-D-D, Filled with water, with outlet on atmosphere \n')
plt.xlabel(r'Time (s)')
plt.legend(['Pipe n°1 $\dot{m}_0$','Pipe n°2 $\dot{m}_1$','Pipe n°3 $\dot{m}_2$','Pipe n°4 $\dot{m}_3$'], fontsize=7,loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
plt.ylabel(r'Mass flow (kg/s)')


# =============================================================================
# Calcul FFT : frequencies analysis and maxima
# =============================================================================

from numpy.fft import fft, fftfreq
import matplotlib.pyplot as plt
import numpy as np


# Your data and FFT calculations
x = list_pressure_1 + 1000 * [0]
N = len(x)
X = fft(x)
freq = fftfreq(N, dt)
N = len(x)
X_abs = np.abs(X[:N // 2]) * 2.0 / N
freq_pos = freq[:N // 2]

# Find the indices of the maximum values above 50 Hz
threshold_freq = 50
max_indices = np.where((freq_pos > threshold_freq) & (X_abs > 0.1 * np.max(X_abs)))[0]
max_frequencies = freq_pos[max_indices]
max_amplitudes = X_abs[max_indices]

plt.figure() # Plot the FFT result
plt.semilogy(freq_pos, X_abs, label="FFT Amplitude")
# Scatter plot of maximum points
plt.scatter(max_frequencies, max_amplitudes, color='red', label='Max Peaks')
# Annotate maximum points with their frequencies
for freq, amp in zip(max_frequencies, max_amplitudes):
    plt.annotate(f'{freq:.2f} Hz', xy=(freq, amp), textcoords='offset points', xytext=(5,5), ha='center', fontsize=8, color='blue')
# Adjust plot properties
plt.xlim(0, 300)
plt.ylim(1000, 1000000)
plt.grid()
plt.xlabel(r"Frequencies (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Fast Fourier Transformation of Pressure in Pipe n°1 $P_1$")
plt.legend()
plt.show()




