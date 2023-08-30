# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:35:32 2023

@author: gauth
"""
# =============================================================================
# README 
#     # idea was to see if convergence with pipe augmentation number, and well it does !
#     so we have a convergence in the pipe
#     the code is slow maybe could be optimised
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

l_pipe = 0.2
l_tank_1 = 10
l_tank_2 = 10

d_tank_1 = 1
d_tank_2 = 1
d_pipe = 0.5

A_pipe = (d_pipe/2)**2 * 3.14
A_tank_1 = (d_tank_1/2)**2 * 3.14
A_tank_2 = (d_tank_2/2)**2 * 3.14

V_tank_1 = l_tank_1 * A_tank_1
V_tank_2 = l_tank_2 * A_tank_2
V_pipe = l_pipe * A_pipe


# =============================================================================
# Physical Constants Input
# =============================================================================
t = 0
dt = 0.0000001

rho0i = 1197
rho1i = 1197
rho2i = 1197

type_simulation = 1 # 0 : gas & 1 : liquid

gamma = 1 # isentropic coefficient
nu = 1e-5 # cinematic viscosity water
eta = 1e-3 # dynamic viscosity water

zeta_geo_1 = 0.5625 #loss coefficient due to geometry
zeta_geo_20 = 9
zeta_geo_i = 0

material_pipe = 50e9 #Pa , Young Modulus aluminium
beta = 2.2e9 #Pa, bulk modulus water


# =============================================================================
# Initialisation 
# =============================================================================

P0i = 0.6e6
P_pipe_i = 0.4e6
P20i = 0.5e6

m0i = 0
m1i = 0
m20i = 0

m0_prev = m0i

m1_prev = m1i
m2_prev = m1i
m3_prev = m1i
m4_prev = m1i
m5_prev = m1i
m6_prev = m1i
m7_prev = m1i
m8_prev = m1i
m9_prev = m1i
m10_prev = m1i
m11_prev = m1i
m12_prev = m1i
m13_prev = m1i
m14_prev = m1i
m15_prev = m1i
m16_prev = m1i
m17_prev = m1i
m18_prev = m1i
m19_prev = m1i

m20_prev = m20i

P0 = P0i + rho0i*10*l_tank_1

P1 = P_pipe_i
P2 = P_pipe_i
P3 = P_pipe_i
P4 = P_pipe_i
P5 = P_pipe_i
P6 = P_pipe_i
P7 = P_pipe_i
P8 = P_pipe_i
P9 = P_pipe_i
P10 = P_pipe_i
P11 = P_pipe_i
P12 = P_pipe_i
P13 = P_pipe_i
P14 = P_pipe_i
P15 = P_pipe_i
P16 = P_pipe_i
P17 = P_pipe_i
P18 = P_pipe_i
P19 = P_pipe_i

P20 = P20i + rho2i*10*l_tank_2


# =============================================================================
# Lists to save results
# =============================================================================
list_pressure_0 = []
list_pressure_10 = []
list_pressure_20 = []

list_mass_flow_0 = []
list_mass_flow_10 = []
list_mass_flow_20 = []

time = []

# =============================================================================
# Functions
# =============================================================================

def a_sound(Pressure_current,Pressure_init,rho_init):
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        return ( beta / rho_init ) / (1+(d_pipe/0.02)*(beta/material_pipe))      
    
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
# Loop Main 
# =============================================================================


while t<0.05:
    t = t + dt


    # boundary input
    
    m0_prev = m1_prev /( 1 + ( (V_pipe/a_sound(P1, P_pipe_i, rho1i))/(9.81/A_tank_1) ) )
    m0_next = m0_prev
    
    # core input 
    
    P0n = P0 + dt * ( (m0_prev-m1_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P1n = P1 + dt * ( (m1_prev-m2_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P2n = P2 + dt * ( (m2_prev-m3_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P3n = P3 + dt * ( (m3_prev-m4_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P4n = P4 + dt * ( (m4_prev-m5_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P5n = P5 + dt * ( (m5_prev-m6_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P6n = P6 + dt * ( (m6_prev-m7_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P7n = P7 + dt * ( (m7_prev-m8_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P8n = P8 + dt * ( (m8_prev-m9_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P9n = P9 + dt * ( (m9_prev-m10_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P10n = P10 + dt * ( (m10_prev-m11_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P11n = P11 + dt * ( (m11_prev-m12_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P12n = P12 + dt * ( (m12_prev-m13_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P13n = P13 + dt * ( (m13_prev-m14_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P14n = P14 + dt * ( (m14_prev-m15_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P15n = P15 + dt * ( (m15_prev-m16_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P16n = P16 + dt * ( (m16_prev-m17_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P17n = P17 + dt * ( (m17_prev-m18_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P18n = P18 + dt * ( (m18_prev-m19_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )
    P19n = P19 + dt * ( (m19_prev-m20_prev) / (V_pipe ) * a_sound(P1,P_pipe_i,rho1i) )

    # boundary input
    
    P20n = P20 + dt * ( ( (m20_prev) * 9.81 / A_tank_2)  )

    zeta_fluid_1 = 64/Reynolds(m1_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi1 = 0.5 * A_pipe**(-2) * ( zeta_fluid_1 + zeta_geo_1 )
    m1_next = m1_prev + dt * ( (A_pipe/l_pipe) * ( P0 - P1 - ( xi1 / rho(P1,P_pipe_i,rho1i) ) * m1_prev * np.abs(m1_prev) ) ) 
    
    zeta_fluid_20 = 64/Reynolds(m1_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi20 = 0.5 * A_pipe**(-2) * ( zeta_fluid_20 + zeta_geo_20 )
    m20_next = m20_prev + dt * ( (A_pipe/l_pipe) * ( P19 - P20 - ( xi20 / rho(P19,P_pipe_i,rho1i) ) * m20_prev * np.abs(m20_prev) ) ) 
    
    zeta_fluid_i = 64/Reynolds(m2_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m2_next = m2_prev + dt * ( (A_pipe/l_pipe) * ( P1 - P2 - ( xi / rho(P1,P_pipe_i,rho1i) ) * m2_prev * np.abs(m2_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m3_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m3_next = m3_prev + dt * ( (A_pipe/l_pipe) * ( P2 - P3 - ( xi / rho(P2,P_pipe_i,rho1i) ) * m3_prev * np.abs(m3_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m4_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m4_next = m4_prev + dt * ( (A_pipe/l_pipe) * ( P3 - P4 - ( xi / rho(P3,P_pipe_i,rho1i) ) * m4_prev * np.abs(m4_prev) ) )
    zeta_fluid_i = 64/Reynolds(m5_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m5_next = m5_prev + dt * ( (A_pipe/l_pipe) * ( P4 - P5 - ( xi / rho(P4,P_pipe_i,rho1i) ) * m5_prev * np.abs(m5_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m6_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m6_next = m6_prev + dt * ( (A_pipe/l_pipe) * ( P5 - P6 - ( xi / rho(P5,P_pipe_i,rho1i) ) * m6_prev * np.abs(m6_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m7_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m7_next = m7_prev + dt * ( (A_pipe/l_pipe) * ( P6 - P7 - ( xi / rho(P6,P_pipe_i,rho1i) ) * m7_prev * np.abs(m7_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m8_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m8_next = m8_prev + dt * ( (A_pipe/l_pipe) * ( P7 - P8 - ( xi / rho(P7,P_pipe_i,rho1i) ) * m8_prev * np.abs(m8_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m9_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m9_next = m9_prev + dt * ( (A_pipe/l_pipe) * ( P8 - P9 - ( xi / rho(P8,P_pipe_i,rho1i) ) * m9_prev * np.abs(m9_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m10_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m10_next = m10_prev + dt * ( (A_pipe/l_pipe) * ( P9 - P10 - ( xi / rho(P9,P_pipe_i,rho1i) ) * m10_prev * np.abs(m10_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m11_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m11_next = m11_prev + dt * ( (A_pipe/l_pipe) * ( P10 - P11 - ( xi / rho(P10,P_pipe_i,rho1i) ) * m11_prev * np.abs(m11_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m12_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m12_next = m12_prev + dt * ( (A_pipe/l_pipe) * ( P11 - P12 - ( xi / rho(P11,P_pipe_i,rho1i) ) * m12_prev * np.abs(m12_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m13_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m13_next = m13_prev + dt * ( (A_pipe/l_pipe) * ( P12 - P13 - ( xi / rho(P12,P_pipe_i,rho1i) ) * m13_prev * np.abs(m13_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m14_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m14_next = m14_prev + dt * ( (A_pipe/l_pipe) * ( P13 - P14 - ( xi / rho(P13,P_pipe_i,rho1i) ) * m14_prev * np.abs(m14_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m15_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m15_next = m15_prev + dt * ( (A_pipe/l_pipe) * ( P14 - P15 - ( xi / rho(P14,P_pipe_i,rho1i) ) * m15_prev * np.abs(m15_prev) ) )
    zeta_fluid_i = 64/Reynolds(m16_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m16_next = m16_prev + dt * ( (A_pipe/l_pipe) * ( P15 - P16 - ( xi / rho(P15,P_pipe_i,rho1i) ) * m16_prev * np.abs(m16_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m17_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m17_next = m17_prev + dt * ( (A_pipe/l_pipe) * ( P16 - P17 - ( xi / rho(P16,P_pipe_i,rho1i) ) * m17_prev * np.abs(m17_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m18_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m18_next = m18_prev + dt * ( (A_pipe/l_pipe) * ( P17 - P18 - ( xi / rho(P17,P_pipe_i,rho1i) ) * m18_prev * np.abs(m18_prev) ) ) 
    zeta_fluid_i = 64/Reynolds(m19_prev,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi = 0.5 * A_pipe**(-2) * ( zeta_fluid_i + zeta_geo_i )
    m19_next = m19_prev + dt * ( (A_pipe/l_pipe) * ( P18 - P19 - ( xi / rho(P18,P_pipe_i,rho1i) ) * m19_prev * np.abs(m19_prev) ) ) 

    # =========================================================================
    #     Loop ton integrate
    # =========================================================================
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
    P10 = P10n
    P11 = P11n
    P12 = P12n
    P13 = P13n
    P14 = P14n
    P15 = P15n
    P16 = P16n
    P17 = P17n
    P18 = P18n
    P19 = P19n
    P20 = P20n

    m0_prev = m0_next
    m1_prev = m1_next
    m2_prev = m2_next
    m3_prev = m3_next
    m4_prev = m4_next
    m5_prev = m5_next
    m6_prev = m6_next
    m7_prev = m7_next
    m8_prev = m8_next
    m9_prev = m9_next
    m10_prev = m10_next
    m11_prev = m11_next
    m12_prev = m12_next
    m13_prev = m13_next
    m14_prev = m14_next
    m15_prev = m15_next
    m16_prev = m16_next
    m17_prev = m17_next
    m18_prev = m18_next
    m19_prev = m19_next
    m20_prev = m20_next


    # =========================================================================
    #     Loop to save 
    # =========================================================================
    list_pressure_0.append(P0n)
    list_pressure_10.append(P10n)
    list_pressure_20.append(P20n)

    list_mass_flow_0.append(m0_prev)
    list_mass_flow_10.append(m10_prev)
    list_mass_flow_20.append(m20_prev)

    time.append(t)
    
    
    # =========================================================================
    #     Print Values
    # =========================================================================
    
    if np.abs(t-0.1)<0.5*dt:
        aa = np.sqrt(a_sound(P10,P_pipe_i,rho1i))
        
        print(f'time : {t} \n P0 : {P0}\n P10 : {P10} \n P20 : {P20}\n m0 : {m0_prev} \n m10 : {m10_prev}\n m20 : {m20_prev} \n a : {aa}')

        aa = np.sqrt(a_sound(P10,P_pipe_i,rho1i))
        print(f'a speed sound non corrected : {aa}')
        aa=aa*np.sqrt(1+(l_pipe/0.1)*(2.2e9/50e9))
        print(f'a speed sound corrected : {aa}')



# =============================================================================
# Plots 
# =============================================================================

plt.figure()
plt.plot(time[::1000],list_pressure_0[::1000])
plt.plot(time[::1000],list_pressure_10[::1000])
plt.plot(time[::1000],list_pressure_20[::1000])

plt.title(f'Liquid Oxygen model : Tank Pipe Tank \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P0','P1','P2'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time[::1000],list_mass_flow_0[::1000])
plt.plot(time[::1000],list_mass_flow_10[::1000])
plt.plot(time[::1000],list_mass_flow_20[::1000])

plt.title(f'mass evolution in filling of Liquid Oxygen of tanks \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m0','m01','m1','m02','m2'])
plt.ylabel(r'$mass flow$ (kg/s)')

# =============================================================================
# FFT Analysis
# =============================================================================

from numpy.fft import fft, fftfreq

# Calcul FFT
x=list_pressure_10 + 10000*[0]
N= len(x)
X = fft(x)  # Transformée de fourier
freq = fftfreq(N, dt)  # Fréquences de la transformée de Fourier
# Calcul du nombre d'échantillon
N= len(x)
# On prend la valeur absolue de l'amplitude uniquement pour les fréquences positives et normalisation
X_abs = np.abs(X[:N//2])*2.0/N
# On garde uniquement les fréquences positives
freq_pos = freq[:N//2]

plt.figure()
plt.plot(freq_pos, X_abs, label="Amplitude absolue")
plt.xlim(0, 200)  # On réduit la plage des fréquences à la zone utile
plt.grid()
plt.xlabel(r"Fréquence (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Transformée de Fourier")
plt.show()
