# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 13:23:32 2023

@author: gauth
"""

# =============================================================================
# README
#     We tried to implement the model of combustion discussed in the report 
#     we used several improvement, which could be seen in the drafts codes
#   ==> see the implicit Euler to have the better synthesis as this code is not functionnal !!!
# The Mixture Ratio is there implemented  which is the most important and interesting aspect of this code, as detailed in the report
# to be fully functionnal, this code would need to be corrected with propellant cosntants instead of water constants, which can be found in the implicit version
# =============================================================================


# =============================================================================
# Imported Modules 
# =============================================================================

import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from scipy.optimize import fsolve
import time 
import cantera as ct 
import numpy as np


# =============================================================================
# Input 
# =============================================================================


a = 1355**2 # speed of sound, but water while we need propellants !
rho = 1197 # density but water while we need propellants !!!! again cf README for explanations


l_pipe = 2
l_combustor = 0.2

d_pipe = 0.01
d_combustor = 0.15 

A_pipe = (d_pipe/2)**2 * 3.14

V_pipe = l_pipe * A_pipe

V_combustor = (0.15/2)**2 * 3.14 *0.2
R = 8.734
R_cc = R/16
A_throat = 0.01

T3 = 500 #K

t = 0
dt = 1e-5

type_simulation = 1
gamma = 1.3
mu = 1e-5
eta = 1e-3
tau = 0.03

P1i = 1e6
P2i = 1e6
P3i = 0.1e6

m1i = 0.0001
m2i = 0.0001
m3i = 0.0001
moxi = m1i * dt
mfui = m2i * dt
mprodi = 0.

m1 = m1i
m2 = m2i
m3 = m3i

mox = moxi
mfu = mfui
mprod = mprodi

P1 = P1i
P2 = P2i
P3 = P3i

# =============================================================================
# Functions 
# =============================================================================

def a_sound(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        # return ( 2.2e9 / rho_init ) / (1+(l_pipe/0.1)*(2.2e9/50e9))
        return 2.2e9 / rho_init

def calculate_rho(Pressure_current,Pressure_init,rho_init):
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

# =============================================================================
# !!! Important !!!
# # mixture ratio function, which gives us the fraction of gas burned, which is then given to cantera for a temperature estimation
# # problem maybe: cantera does not do full combustion of this ! so maybe does cantera already do this operation and that is why it was not implemented in the implicit version
# # as a matter of fact, cantera do partial combustion and stoechiometric equilibrium 
# =============================================================================

def mu_compos(m1,m2): # model explained in the report in the section about chamber and mixture ratio 
    if m2 == 0 :
        mu_ox_burned, mu_ox_out, mu_fu_burned, mu_fu_out = 1,1,1,1
    else:
        MRin = m1/m2
        # MRst = 0.5 * 32/2
        MRst = 2 
        if MRin <= MRst:
            mu_ox_burned = MRin / (MRin + 1)
            mu_ox_out = 0
            mu_fu_burned = MRin / ( (MRin + 1) * MRst)
            mu_fu_out = (1 - MRin/MRst)/(MRin + 1)
        else :
            mu_ox_burned = MRst/ (MRin + 1)
            mu_ox_out =(MRin - MRst)/(MRin + 1)
            mu_fu_burned = 1 / (1 + MRin)
            mu_fu_out = 0
    return [mu_ox_burned, mu_ox_out, mu_fu_burned, mu_fu_out]


list_pressure_1 = []
list_pressure_2 = []
list_pressure_3 = []

list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []

list_temperature_3 = []

list_mass_tot = []
list_mass_ox = []
list_mass_fu = []
list_mass_prod = []
list_mixture_ratio_inside = []

time_storage = []

start = time.time() 
print("Start computing") 

gas1 = ct.Solution('gri30.yaml')

while t<0.6:
    
    ############################################################  
    #### LOOP AND CALCULUS 
    ############################################################
    
    t = t + dt
#    # useless: was part of the drafts code versions, where we needed to compute the ratio of burned and out species, while we do it now with the mixture ratio function
#    # mox = (m1 * (1 - mu_compos(m1, m2)[0]) - m3*mu_compos(m1, m2)[1] ) * dt
#    # mfu = (m2 * (1 - mu_compos(m1, m2)[2]) - m3*mu_compos(m1, m2)[3] ) * dt 
#    # mprod = ( m1 * mu_compos(m1, m2)[0] +  m2 * mu_compos(m1, m2)[2] - ( 1 - mu_compos(m1, m2)[1] - mu_compos(m1, m2)[3] )*m3 ) * dt

    # we compute energy we release thanks to input without delay and with what is consumed with total combustion 
    # we use perez rocca 
    gas1.Y = f' H2:{ mu_compos(m1, m2)[2] * m2 } , O2:{ mu_compos(m1, m2)[0] * m1 }'
    gas1.TP = 300, P3
    gas1.equilibrate('HP')
    T3n = gas1.T
    P1n = P1
    P2n = P2
    m1n = m1 + dt * A_pipe/l_pipe * (P1 - P3)
    m2n = m2 + dt * A_pipe/l_pipe * (P2 - P3)
    # P3n = P3 + dt * R_cc * T3 / V_combustor * ( m1 + m2 - m3) # when temperature constant
    # and when temperature not constant, we have : 
    P3n = P3 + dt * R_cc * T3 / V_combustor * ( m1 + m2 - m3 + V_combustor * P3 / R_cc * (T3n - T3)/dt * (+1 * 1/T3**2) )
    m3n = A_throat * np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) / np.sqrt(R_cc * T3n) * P3n
    P1 = P1n
    P2 = P2n
    P3 = P3n
    m1 = m1n
    m2 = m2n
    m3 = m3n
    T3 = T3n
    
    ############################################################  
    #### SAVE 
    ############################################################
    
    list_pressure_1.append(P1)
    list_pressure_2.append(P2)
    list_pressure_3.append(P3)

    list_mass_flow_1.append(m1)
    list_mass_flow_2.append(m2)
    list_mass_flow_3.append(m3)
    
    list_temperature_3.append(T3)
    
    time_storage.append(t)
    
    list_mass_tot.append(mox + mfu + mprod)
    list_mass_ox.append(mox)
    list_mass_fu.append(mfu)
    list_mass_prod.append(mprod)
    list_mixture_ratio_inside.append( (mu_compos(m1, m2)[0] * m1) / (mu_compos(m1, m2)[2] * m2) )
    
end = time.time() 
print(f'time of simulation : {end - start}')   

 
# =============================================================================
# # Plots and graphs
# =============================================================================


plt.figure()
plt.plot(time_storage[::2],list_pressure_1[::2])
plt.plot(time_storage[::2],list_pressure_2[::2])
plt.plot(time_storage[::2],list_pressure_3[::2])
plt.title(f'Pressure evolution Combustor fed at changing temperature \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P1 tank supply ox','P2 tank supply fuel','P3 combustion chamber pressure'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time_storage[::2],list_mass_flow_1[::2])
plt.plot(time_storage[::2],list_mass_flow_2[::2])
plt.plot(time_storage[::2],list_mass_flow_3[::2])

plt.title(f'massflow evolution in Combustor fed at changing temperature \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m1 ox','m2 fuel','m3 out'])
plt.ylabel(r'$mass flow$ (kg/s)')
# plt.ylim(-0, 50)

plt.figure()
plt.plot(time_storage[::2],list_temperature_3[::2]) 
plt.title(f'temperature evolution in Combustor fed at changing temperature \n')
plt.xlabel(r'$time$ (s)')
plt.ylabel(r'$T$ (K)')
# plt.ylim(T3 - 1000 , T3 + 1000)

plt.figure()
plt.plot(time_storage[::2],list_mixture_ratio_inside[::2]) 
plt.title(f'mixture ratio incoming evolution in Combustor fed at changing temperature \n')
plt.xlabel(r'$time$ (s)')
plt.ylabel(r'$MR$ (%)')
plt.ylim(-0, 10)


plt.figure()
plt.plot(time_storage[::2],list_mass_tot[::2]) 
plt.plot(time_storage[::2],list_mass_ox[::2]) 
plt.plot(time_storage[::2],list_mass_fu[::2]) 
plt.plot(time_storage[::2],list_mass_prod[::2]) 
plt.title(f'mixture inside evolution in Combustor fed at changing temperature \n')
plt.xlabel(r'$time$ (s)')
plt.ylabel(r'$mass$ (kg)')
plt.legend(['total mass','m1 ox','m2 fuel','m products'])
# plt.ylim(-0, 1)

