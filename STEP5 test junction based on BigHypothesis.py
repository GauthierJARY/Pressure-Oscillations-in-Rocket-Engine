# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 16:41:12 2023

@author: gauth
"""

"""
Comment :
    test to find a solution to join and connect => that is to say to have a loss between pipe and tank 
    the fact is it is already taken with zeta into account 
    not interesting code 
    this was done when i first discovered lumped parameter hypothesis but is no longer accurate 
"""

import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys

###################
# geometric input 
###################

l_pipe_1 = 2
l_pipe_2 = 2
l_tank_1 = 10
l_tank_2 = 10

d_tank_1 = 1
d_tank_2 = 1
d_pipe_1 = 0.5
d_pipe_2 = 0.5

A_pipe_1 = (d_pipe_1/2)**2 * 3.14
A_pipe_2 = (d_pipe_2/2)**2 * 3.14
A_tank_1 = (d_tank_1/2)**2 * 3.14
A_tank_2 = (d_tank_2/2)**2 * 3.14

V_tank_1 = l_tank_1 * A_tank_1
V_tank_2 = l_tank_2 * A_tank_2
V_pipe_1 = l_pipe_1 * A_pipe_1
V_pipe_2 = l_pipe_2 * A_pipe_2

###################
# initial input 
###################
Ppi = 0.6e6 
P0i = 0.4e6
P1i = 0.4e6
P2i = 0.4e6
Pppi = 0.5e6

mpi = 0
m0i = 0
m1i = 0
m2i = 0
m3i = 0
mppi = 0

t = 0
dt = 0.0000005
# dt = 1e-3 converge (vir la différence entre deux division du pas si elle est minime le resultat a converge)
mp_prev = mpi
m0_prev = m0i
m1_prev = m1i
m2_prev = m2i
mpp_prev = mppi

rhopi = 1197
rho0i = 1197
rho1i = 1197
rho2i = 1197
rhoppi = 1197

type_simulation = 1
# 0 : gas & 1 : liquid

Pp = Ppi + rhopi*10*l_tank_1
P0 = P0i
P1 = P1i
P2 = P2i
Ppp = Pppi + rhoppi*10*l_tank_2

###################
# constant input 
###################

gamma = 1
mu = 1e-5
zeta1 = 0
zeta2 = 0
K1 = 0.16
K2 = 0.5

list_Pressure_p = []
list_Pressure_0 = []
list_Pressure_1 = []
list_Pressure_2 = []
list_Pressure_pp = []

list_mass_flow_p = []
list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_pp = []

time = []


def a_sound(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        # return ( 2.2e9 / rho_init ) / (1+(l_pipe_1/0.1)*(2.2e9/50e9))
        return ( 2.2e9 / rho_init )

def rho(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 : 
        return rho_init * (Pressure_current/Pressure_init)**(1/gamma)
    else :
        return rho_init

def Reynolds(m_flow, length):
    m=m_flow
    if m_flow==0:
        m=1e-12
    return m / (length*mu)
    # return 300

# while np.abs(P0-P2)>5e-5 :
while t<0.012:
    t = t + dt
    # sys.stdout.write("\rSimulation time advancement: {0:.2f} %".format(100* t/3))
    # sys.stdout.flush()
    # loop to calculus
     
    Ppn = Pp + dt * ( 9.81 / A_tank_1 ) * (-mp_prev)
    m0_prev = K1 *np.sign(P0-Pp)*np.sqrt(2*rho(P0,P0i,rho0i)*np.abs(P0-Pp))
    P0n = P0 + dt * ( (m0_prev-m1_prev) / (V_pipe_1) * a_sound(P0,P0i,rho0i) )    
    P1n = P1 + dt * ( (m1_prev-m2_prev) / (V_pipe_1) * a_sound(P1,P1i,rho1i) )
    # mpp_prev = K2 * np.sqrt( 2*rho(P2,P2i,rho2i)*(Ppp-P2) )
    Pppn = Ppp + dt * ( (  mpp_prev* 9.81 / A_tank_2)  )

    xi1 = 64/Reynolds(m1_prev,l_pipe_1)*1/d_pipe_1 * 0.5 * A_pipe_1**(-2) + zeta1 * 0.5 * A_pipe_1**(-2)
    m1_next = m1_prev + dt * ( (A_pipe_1/l_pipe_1) * ( P0 - P1 - ( xi1 / rho(P1,P1i,rho1i) ) * m1_prev * np.abs(m1_prev) ) ) 
    
    xi2 = 64/Reynolds(m2_prev,l_pipe_2)*1/d_pipe_2 * 0.5 * A_pipe_2**(-2) + zeta2 * 0.5 * A_pipe_2**(-2)
    m2_next = m2_prev + dt * ( (A_pipe_1/l_pipe_2) * ( P1 - P2 - ( xi2 / rho(P1,P1i,rho1i) ) * m2_prev * np.abs(m2_prev) ) ) 
    P2n = -( (m2_next/K2)**2 / (2*rho(P2,P2i,rho2i))-Ppp)
    
    # loop to integrate 
    Ppn = Pp
    Pppn = Ppp
    P0 = P0n
    P1 = P1n
    P2 = P2n
    m1_prev = m1_next
    m2_prev = m2_next
    mp_prev = mp_prev # 0D element
    mpp_prev = mpp_prev # 0D element

    
    # loop to save 
    list_Pressure_p.append(Ppn)
    list_Pressure_0.append(P0n)
    list_Pressure_1.append(P1n)
    list_Pressure_2.append(P2n)
    list_Pressure_pp.append(Pppn)
    list_mass_flow_0.append(m0_prev)
    list_mass_flow_1.append(m1_prev)
    list_mass_flow_2.append(m2_prev)
    list_mass_flow_pp.append(mp_prev)
    list_mass_flow_p.append(mpp_prev)

    time.append(t)
    
    if np.abs(t-0.1)<dt:
        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        
        print(f'time : {t} \n P0 : {P0}\n P1 : {P1} \n P2 : {P2}\n m0 : {m0_prev} \n m1 : {m1_prev}\n m2 : {mp_prev} \n a : {aa}')

        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        print(f'a speed sound non corrected : {aa*np.sqrt(1+(l_pipe_1/0.1)*(2.2e9/50e9))}')
        print(f'a speed sound corrected : {aa}')
        print(f'harmonic frequencies : {aa/(2*2*l_pipe_1)}')


plt.figure()
plt.plot(time[::100],list_Pressure_0[::100])
plt.plot(time[::100],list_Pressure_p[::100])
plt.plot(time[::100],list_Pressure_1[::100])
plt.plot(time[::100],list_Pressure_pp[::100])
plt.plot(time[::100],list_Pressure_2[::100])
plt.title(f'Liquid Oxygen model : Tank Pipe Tank \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P0','Pp','P1', 'Ppp','P2'])
plt.ylabel(r'$Pressure$ (Pa)')

# plt.figure()
# plt.plot(time[::1000],Q[::1000])
# plt.plot(time[::1000],Qr[::1000])
# plt.plot(time[::1000],Qre[::1000])
# plt.plot(time[::1000],Qrez[::1000])

# plt.title(f'mass evolution in filling of Liquid Oxygen of tanks \n')
# plt.xlabel(r'$time$ (s)')
# plt.legend(['m0','m1_prev','m2','m3'])
# plt.ylabel(r'$mass flow$ (kg/s)')



from numpy.fft import fft, fftfreq
plt.figure()
# Calcul FFT
x=list_Pressure_1
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
plt.grid()
plt.xlabel(r"Fréquence (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Transformée de Fourier")
plt.show()

# Get the indices of maximum element in numpy array
result = np.where(X_abs == np.amax(X_abs))

