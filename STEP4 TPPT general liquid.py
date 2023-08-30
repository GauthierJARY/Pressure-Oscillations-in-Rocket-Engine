# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:14:36 2023

@author: gauth
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

P0i = 0.6e6
P1i = 0.4e6
P2i = 0.5e6

m0i = 0
m1i = 0
m2i = 0
m3i = 0

t = 0
dt = 0.00000001

m0_prev = m0i
m1_prev = m1i
m2 = m2i
m3_prev = m3i


rho0i = 1197
rho1i = 1197
rho2i = 1197

type_simulation = 1
# 0 : gas & 1 : liquid

P0 = P0i + rho0i*10*l_tank_1
P1 = P1i
P2 = P2i + rho2i*10*l_tank_2

###################
# constant input 
###################

gamma = 1
mu = 1e-5
zeta1 = 0.5625
zeta2 = 9

P = []
Pr = []
Pre = []

Q = []
Qr = []
Qre = []
Qrez = []

mass_tot = []

time = []


def a_sound(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        return 2.2e9 / rho_init

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
while t<0.5:
    t = t + dt
    # sys.stdout.write("\rSimulation time advancement: {0:.2f} %".format(100* t/3))
    # sys.stdout.flush()
    # loop to calculus
     
    P0n = P0 + dt * ( a_sound(P0,P0i,rho0i) / V_tank_1 ) * (m0_prev - m1_prev)
    P1n = P1 + dt * ( (m1_prev-m2) / (V_pipe_1 + V_pipe_2) * a_sound(P1,P1i,rho1i) )
    P2n = P2 + dt * ( a_sound(P2,P2i,rho2i) / V_tank_2 ) * (m2 - m3_prev)

    xi1 = 64/Reynolds(m1_prev,l_pipe_1)*1/d_pipe_1 * 0.5 * A_pipe_1**(-2) + zeta1 * 0.5 * A_pipe_1**(-2)
    m1_next = m1_prev + dt * ( 1/l_pipe_1 * ( P0*A_pipe_1 - P1 * A_pipe_1 - xi1 / rho(P1,P1i,rho1i) * A_pipe_1 * m1_prev * np.abs(m1_prev) ) ) 
    
    xi2 = 64/Reynolds(m2,l_pipe_2)*1/d_pipe_2 * 0.5 * A_pipe_2**(-2) + zeta2 * 0.5 * A_pipe_2**(-2)
    m2n = m2 + dt * ( 1/l_pipe_2 * ( P1*A_pipe_1 - P2 * A_pipe_1 - xi2 / rho(P1,P1i,rho1i) * A_pipe_2 * m2 * np.abs(m2) ) ) 

    m0_next = m0_prev + dt * ( 1/l_tank_1 * (-A_pipe_1*P0) - rho(P0,P0i,rho0i) * 10 * A_tank_1)    
    
    m3_next = m3_prev + dt * ( 1/l_tank_2 * (A_pipe_2*P2) - rho(P2,P2i,rho2i) * 10 * A_tank_2)    

    # loop to integrate 
    P0 = P0n
    P1 = P1n
    P2 = P2n
    m1_prev = m1_next
    m2 = m2n
    # m0_prev = m0_next
    # m3_prev = m3_next
    
    # loop to save 
    P.append(P0n)
    Pr.append(P1n)
    Pre.append(P2n)
    Q.append(m0_prev)
    Qr.append(m1_prev)
    Qre.append(m2)
    Qrez.append(m3_prev)
    time.append(t)
    mass_tot.append( rho(P1,P1i,rho1i)*(V_pipe_1+V_pipe_2) + rho(P0,P0i,rho0i)*(V_tank_1) + rho(P2,P2i,rho2i)*(V_tank_2) )
    
plt.figure()
plt.plot(time,P)
plt.plot(time,Pr)
plt.plot(time,Pre)
plt.title(f'Liquid Oxygen model : Tank Pipe Tank \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P0','P1', 'P2'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time,Q)
plt.plot(time,Qr)
plt.plot(time,Qre)
plt.plot(time,Qrez)

plt.title(f'mass evolution in filling of Liquid Oxygen of tanks \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m0','m1_prev','m2','m3'])
plt.ylabel(r'$mass flow$ (kg/s)')



from numpy.fft import fft, fftfreq
plt.figure()
# Calcul FFT
x=Pr
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

print('Returned tuple of arrays :', result)
print('List of Indices of maximum element :', result[0])
