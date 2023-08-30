# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 13:00:36 2023

@author: gauth
"""

""" 
Comment : 
    made to compare implicit and explicit with STEP4 A1 model of Chiara 
    Lumped Hyopthesis
    nothing new, should be identical to Chiara s results
    too much dampling it seems compared to what was observed 
    dt can be adjusted as it is implicit so more robust to divergence and iterations 
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

t = 0
dt = 0.000001
# dt = 1e-3 converge (vir la différence entre deux division du pas si elle est minime le resultat a converge)

m0_prev = m0i
m1_prev = m1i
m2_prev = m2i

rho0i = 1197
rho1i = 1197
rho2i = 1197

P0 = P0i + rho0i*10*l_tank_1
P1 = P1i
P2 = P2i + rho2i*10*l_tank_2

type_simulation = 1
# 0 : gas & 1 : liquid


###################
# constant input 
###################

gamma = 1
mu = 1e-5
eta = 1e-3
zeta1 = 0.5625
zeta2 = 9

list_pressure_0 = []
list_pressure_1 = []
list_pressure_2 = []

list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []

list_energy = []
list_power = []
list_tot =[]

time = []


def a_sound(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        # return ( 2.2e9 / rho_init ) / (1+(l_pipe_1/0.1)*(2.2e9/50e9))
        return 2.2e9 / rho_init

def rho_calculate(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 : 
        return rho_init * (Pressure_current/Pressure_init)**(1/gamma)
    else :
        return rho_init

def Reynolds(m_flow, length, section):
    m=m_flow
    if m_flow==0:
        m=1e-12
    # return m / (length*mu)
    return (m * length ) / ( section * eta )

a= 1355**2
rho = 1355
from scipy.optimize import fsolve

def function(S):
    C = a/V_pipe_1
    Cprime = A_tank_1/9.81
    C=1/C
    Cprime = 1/Cprime
    m0 = 1/(1+C/Cprime)*m1
    lambda1 = 64/Reynolds(m1,l_pipe_1,A_pipe_1)*1/d_pipe_1 *l_pipe_1
    xi1 = 0.5 * A_pipe_1**(-2) * ( lambda1 + zeta1 )
    lambda2 = 64/Reynolds(m2,l_pipe_1,A_pipe_1)*1/d_pipe_1 *l_pipe_1
    xi2 = 0.5 * A_pipe_1**(-2) * ( lambda2 + zeta2 )
    return np.array([
        S[0] - P0 - dt * a / V_pipe_1 *(m0 - S[3]) ,
        S[1] - P1 - dt * a / V_pipe_2 *(S[3] - S[4]) ,
        S[2] - P2 - dt * 9.81/A_tank_2 * S[4] ,
        S[3] - m1 - dt * ( (A_pipe_1/l_pipe_1) * ( S[0] - S[1] - ( xi1 / rho ) * S[3] * np.abs(S[3]) ) ) ,
        S[4] - m2 - dt * ( (A_pipe_2/l_pipe_2) * ( S[1] - S[2] - ( xi2 / rho ) * S[4] * np.abs(S[4]) ) ) 
        ])




# while np.abs(P0-P2)>5e-5 :
while t<0.2:
    t = t + dt
    root = fsolve(function, [1, 1,1,1,1])
    P0 = root[0]
    P1 = root[1]
    P2 = root[2]
    m1 = root[3]
    m2 = root[4]

    C = a/V_pipe_1
    Cprime = A_tank_1/9.81
    C=1/C
    Cprime = 1/Cprime
    m0 = 1/(1+C/Cprime)*m1

    # loop to save 
    list_pressure_0.append(P0)
    list_pressure_1.append(P1)
    list_pressure_2.append(P2)

    list_mass_flow_0.append(m0)
    list_mass_flow_1.append(m1)
    list_mass_flow_2.append(m2)
    
    time.append(t)

    
    if np.abs(t-0.1)<0.5*dt:
        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        
        print(f'time : {t} \n P0 : {P0}\n P1 : {P1} \n P2 : {P2}\n m0 : {m0_prev} \n m1 : {m1_prev}\n m2 : {m2_prev} \n a : {aa}')

        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        print(f'a speed sound non corrected : {aa}')
        # aa=aa/np.sqrt(1+(l_pipe_1/0.1)*(2.2e9/50e9))
        print(f'a speed sound corrected : {aa}')
        V1=V_pipe_1
        lambda1 = 64/Reynolds(m1,l_pipe_1,A_pipe_1)*1/d_pipe_1 *l_pipe_1
        xi1 = 0.5 * A_pipe_1**(-2) * ( lambda1 + zeta1 )
        lambda2 = 64/Reynolds(m2,l_pipe_1,A_pipe_1)*1/d_pipe_1 *l_pipe_1
        xi2 = 0.5 * A_pipe_1**(-2) * ( lambda2 + zeta2 )
        lin1 = A_pipe_1/l_pipe_1*xi1/rho*np.abs(m1)
        # lin1 = A_pipe_1/l_pipe_1*xi1/rho1i * 1/1e-6
        lin2 = A_pipe_2/l_pipe_2*xi2/rho*np.abs(m2)
        # lin2 = A_pipe_2/l_pipe_2*xi2/rho2i * 1/1e-6

        A = np.array( [ [0,0,0,-9.81/A_tank_1,0],[0,0,0,aa**2/V1,-aa**2/V1],[0,0,0,0,9.81/A_tank_2],[A_pipe_1/l_pipe_1,-A_pipe_1/l_pipe_1,0,lin1,0],[0,A_pipe_2/l_pipe_2,-A_pipe_2/l_pipe_2,0,lin2]])
        vap = np.linalg.eigvals(A)
        print(f'eigen values are : {vap}')
        freqq = np.imag(vap[0])/(2*3.14)
        print(f'eigen frequencies are : {freqq} Hz')


plt.figure()
plt.plot(time[::1000],list_pressure_0[::1000])
plt.plot(time[::1000],list_pressure_1[::1000])
plt.plot(time[::1000],list_pressure_2[::1000])

plt.title(f'Liquid Oxygen model : Tank Pipe Tank \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P0','P1','P2'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time[::1000],list_mass_flow_0[::1000])
plt.plot(time[::1000],list_mass_flow_1[::1000])
plt.plot(time[::1000],list_mass_flow_2[::1000])

plt.title(f'mass evolution in filling of Liquid Oxygen of tanks \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m0','m01','m1','m02','m2'])
plt.ylabel(r'$mass flow$ (kg/s)')



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
plt.grid()
plt.xlabel(r"Fréquence (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Transformée de Fourier")
plt.show()

# Get the indices of maximum element in numpy array
result = np.where(X_abs == np.amax(X_abs))

