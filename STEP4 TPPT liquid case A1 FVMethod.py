# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:01:53 2023

@author: gauth
"""
"""
Comment : 
    
    # not an interesting code because of line 140 where P1 = P_next + P_before / 2
    that was a try on FV and try to approximate the value of the boundary but not accurate !! (at least with so few elements)
    not a usefull code for the report
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
dt = 0.000001
# dt = 1e-3 converge (vir la différence entre deux division du pas si elle est minime le resultat a converge)

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

P01 = P1i
P02 = P1i 
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
        return ( 2.2e9 / rho_init ) / (1+(l_pipe_1/0.1)*(2.2e9/50e9))
        # return 2.2e9 / rho_init

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
while t<0.2:
    t = t + dt
    # sys.stdout.write("\rSimulation time advancement: {0:.2f} %".format(100* t/3))
    # sys.stdout.flush()
    # loop to calculus
    
    P01n = P01 + dt * ( - a_sound(P01, P1i, rho1i) / V_pipe_1 * ((m2)/2) )
    P02n = P02 + dt * ( - a_sound(P02, P1i, rho1i) / V_pipe_2 * (-(m1_prev)/2) )
    

    P0n = P0 + dt * ( 9.81 / A_tank_1 ) * (m0_prev - m1_prev)
    # P1n = P1 + dt * ( (m1_prev-m2) / (V_pipe_1 ) * a_sound(P1,P1i,rho1i) )
    P2n = P2 + dt * ( ( (m2 - m3_prev) * 9.81 / A_tank_2)  )

    xi1 = 64/Reynolds(m1_prev,l_pipe_1)*1/d_pipe_1 * 0.5 * A_pipe_1**(-2) + zeta1 * 0.5 * A_pipe_1**(-2)
    P1 = (P01 + P02) / 2
    m1_next = m1_prev + dt * ( (A_pipe_1/l_pipe_1) * ( P0 - P1 - ( xi1 / rho(P1,P1i,rho1i) ) * m1_prev * np.abs(m1_prev) ) ) 
    
    xi2 = 64/Reynolds(m2,l_pipe_2)*1/d_pipe_2 * 0.5 * A_pipe_2**(-2) + zeta2 * 0.5 * A_pipe_2**(-2)
    m2n = m2 + dt * ( (A_pipe_1/l_pipe_2) * ( P1 - P2 - ( xi2 / rho(P1,P1i,rho1i) ) * m2 * np.abs(m2) ) ) 

    m0_next = m0_prev + dt * ( 1/l_tank_1 * (-A_pipe_1*P0) - rho(P0,P0i,rho0i) * 10 * A_tank_1)    
    
    m3_next = m3_prev + dt * ( 1/l_tank_2 * (A_pipe_2*P2) - rho(P2,P2i,rho2i) * 10 * A_tank_2)    

    # loop to integrate 
    P0 = P0n
    # P1 = P1n
    P2 = P2n
    m1_prev = m1_next
    m2 = m2n
    P01 = P01n
    P02 = P02n
    # m0_prev = m0_next
    # m3_prev = m3_next
    
    # loop to save 
    P.append(P0n)
    Pr.append(P1)
    Pre.append(P2n)
    Q.append(m0_prev)
    Qr.append(m1_prev)
    Qre.append(m2)
    Qrez.append(m3_prev)
    time.append(t)
    mass_tot.append( rho(P1,P1i,rho1i)*(V_pipe_1+V_pipe_2) + rho(P0,P0i,rho0i)*(V_tank_1) + rho(P2,P2i,rho2i)*(V_tank_2) )
    if np.abs(t-0.1)<dt:
        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        
        print(f'time : {t} \n P0 : {P0}\n P1 : {P1} \n P2 : {P2}\n m0 : {m0_prev} \n m1 : {m1_prev}\n m2 : {m2} \n a : {aa}')

        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        print(f'a speed sound non corrected : {aa}')
        # aa=aa/np.sqrt(1+(l_pipe_1/0.1)*(2.2e9/50e9))
        print(f'a speed sound corrected : {aa}')
        V1=V_pipe_1*4
        lin1 = A_pipe_1/l_pipe_1*xi1/rho1i*np.abs(m1_prev)
        # lin1 = A_pipe_1/l_pipe_1*xi1/rho1i * 1/1e-6
        lin2 = A_pipe_2/l_pipe_2*xi2/rho2i*np.abs(m2)
        # lin2 = A_pipe_2/l_pipe_2*xi2/rho2i * 1/1e-6

        A = np.array( [ [0,0,0,-9.81/A_tank_1,0],[0,0,0,aa**2/V1,-aa**2/V1],[0,0,0,0,9.81/A_tank_2],[A_pipe_1/l_pipe_1,-A_pipe_1/l_pipe_1,0,lin1,0],[0,A_pipe_2/l_pipe_2,-A_pipe_2/l_pipe_2,0,lin2]])
        vap = np.linalg.eigvals(A)
        print(f'eigen values are : {vap}')
        freqq = np.imag(vap[0])/(2*3.14)
        print(f'eigen frequencies are : {freqq} Hz')


plt.figure()
plt.plot(time[::1000],P[::1000])
plt.plot(time[::1000],Pr[::1000])
plt.plot(time[::1000],Pre[::1000])
plt.title(f'Liquid Oxygen model : Tank Pipe Tank \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['Tank n°1 pressure P0','Pipe pressure P1', ' Tank n°2 pressure P2'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time[::1000],Q[::1000])
plt.plot(time[::1000],Qr[::1000])
plt.plot(time[::1000],Qre[::1000])
plt.plot(time[::1000],Qrez[::1000])

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

