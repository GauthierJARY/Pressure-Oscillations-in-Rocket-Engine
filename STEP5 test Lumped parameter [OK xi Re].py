# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:32:54 2023

@author: gauth
"""
"""
Comment : 
    very IMPORTANT code to me !!
    it is the 1st try once i discovered the lumped hyptohesis 
    try to implement energy loss
    has eigen values etc 
    very close to STEP4 A1 TPPT Chiara 
    i tried to implement new zeta in it 
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
nu = 1e-5
zeta1 = 0.5625
# zeta1 = 1
zeta2 = 9
# zeta2 = 0.5
eta = 1e-3

list_pressure_0 = []
list_pressure_1 = []
list_pressure_2 = []

list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []

time = []


def a_sound(Pressure_current,Pressure_init,rho_init):
    return ( 2.2e9 / rho_init ) / (1+(d_pipe_1/0.02)*(2.2e9/50e9))
    # return 2.2e9 / rho_init

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

# while np.abs(P0-P2)>5e-5 :
while t<0.05:
    t = t + dt

    # loop to calculus
    m0_prev = m1_prev /( 1 + ( (V_pipe_1/a_sound(P1, P1i, rho1i))/(9.81/A_tank_1) ) )
    m0_next = m0_prev
    # P0n = P0 + dt * ( 9.81 / A_tank_1 ) * (- m0_prev)
    P0n = P0 + dt * ( (m0_prev-m1_prev) / (V_pipe_1 ) * a_sound(P1,P1i,rho1i) )
    P1n = P1 + dt * ( (m1_prev-m2_prev) / (V_pipe_1 ) * a_sound(P1,P1i,rho1i) )
    P2n = P2 + dt * ( ( (m2_prev) * 9.81 / A_tank_2)  )

    lambda1 = 64/Reynolds(m1_prev,l_pipe_1,A_pipe_1)*1/d_pipe_1 *l_pipe_1
    xi1 = 0.5 * A_pipe_1**(-2) * ( lambda1 + zeta1 )
    m1_next = m1_prev + dt * ( (A_pipe_1/l_pipe_1) * ( P0 - P1 - ( xi1 / rho(P1,P1i,rho1i) ) * m1_prev * np.abs(m1_prev) ) ) 
    
    lambda2 = 64/Reynolds(m2_prev,l_pipe_2,A_pipe_2)*2/d_pipe_2 *l_pipe_2
    xi2 = 0.5 * A_pipe_2**(-2) * ( lambda2 + zeta2 )
    m2_next = m2_prev + dt * ( (A_pipe_1/l_pipe_2) * ( P1 - P2 - ( xi2 / rho(P1,P1i,rho1i) ) * m2_prev * np.abs(m2_prev) ) ) 
    

    
    # loop to integrate 
    P0 = P0n
    P1 = P1n
    P2 = P2n

    m0_prev = m0_next
    m1_prev = m1_next
    m2_prev = m2_next

    # loop to save 
    list_pressure_0.append(P0n)
    list_pressure_1.append(P1)
    list_pressure_2.append(P2n)

    list_mass_flow_0.append(m0_prev)
    list_mass_flow_1.append(m1_prev)
    list_mass_flow_2.append(m2_prev)
    
    time.append(t)

    
    if np.abs(t-0.1)<0.5*dt:
        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        
        print(f'time : {t} \n P0 : {P0}\n P1 : {P1} \n P2 : {P2}\n m0 : {m0_prev} \n m1 : {m1_prev}\n m2 : {m2_prev} \n a : {aa}')

        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        print(f'a speed sound non corrected : {aa*np.sqrt(1+(l_pipe_1/0.1)*(2.2e9/50e9))}')
        print(f'a speed sound corrected : {aa}')
        V1=V_pipe_1
        lin1 = A_pipe_1/l_pipe_1*xi1/rho1i*np.abs(m1_prev)
        lin2 = A_pipe_2/l_pipe_2*xi2/rho2i*np.abs(m2_prev)        
        # lin1 = A_pipe_1/l_pipe_1*xi1/rho1i * 1/1e-6
        # lin2 = A_pipe_2/l_pipe_2*xi2/rho2i * 1/1e-6

        A = np.array( [ [0,0,0,-9.81/A_tank_1,0],[0,0,0,aa**2/V1,-aa**2/V1],[0,0,0,0,9.81/A_tank_2],[A_pipe_1/l_pipe_1,-A_pipe_1/l_pipe_1,0,lin1,0],[0,A_pipe_2/l_pipe_2,-A_pipe_2/l_pipe_2,0,lin2]])
        vap = np.linalg.eigvals(-A)
        print(f'eigen values are : {vap}')
        freqq = np.imag(vap[0])/(2*3.14)
        print(f'eigen frequencies are : {freqq} Hz')
        print(f'harmonic frequencies : {aa/(2*2*l_pipe_1)}')



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



from numpy.fft import fft, fftfreq
plt.figure()
# Calcul FFT
x=list_pressure_1
x=x+ [0]*int(len(x)/8)
N= len(x)
X = fft(x)  # Transformée de fourier
freq = fftfreq(N, dt)  # Fréquences de la transformée de Fourier
# Calcul du nombre d'échantillon
# On prend la valeur absolue de l'amplitude uniquement pour les fréquences positives et normalisation
X_abs = np.abs(X[:N//2])*2.0/N
# On garde uniquement les fréquences positives
freq_pos = freq[:N//2]
m=max(X_abs)
# index=X_abs.index(m)
plt.plot(freq_pos, X_abs, label="Amplitude absolue")
plt.ylim(0,0.2e6)
plt.xlim(0, 200)  # On réduit la plage des fréquences à la zone utile
plt.grid()
plt.xlabel(r"Frequency (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Fast Fourier Transformation of the pressure $P_1$")
plt.show()

# Get the indices of maximum element in numpy array
result = np.where(X_abs == np.amax(X_abs))

