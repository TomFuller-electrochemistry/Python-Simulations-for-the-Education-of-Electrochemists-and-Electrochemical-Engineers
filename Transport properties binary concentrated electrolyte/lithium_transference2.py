# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:40:57 2020

@author: tomff
"""
import numpy as np
from scipy.constants import e, R, N_A
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#from matplotlib.widgets import Slider

L = 20 #thickness of separator, microns
cd = 350 #current density, A/m^2
F = e*N_A # Faraday's constant, C/eq
c0 = 1000 #concentration of electrolyte, mol/m^3
eps = 0.3
t_plus = 0.0 #lithium ion transference number
D_b = 2.79e-10 #diffusivity, m^2/s
D_e = eps**1.5 * D_b


def conductivity(c):
    k_eff = 0.0911 + 1.9101e-3*c - 1.05e-6*c*c + 0.1554e-9*c*c*c
    return k_eff

def dy_dx(y,x):
    c = ((cd/F)*(1-t_plus)/D_e)*(L/2-x)*1e-6 + c0
    k_eff = eps**1.5 * (0.0911 + 1.9101e-3*c - 1.05e-6*c*c + 0.1554e-9*c*c*c)
    return -1e-6*cd/k_eff

#SET UP PLOT
plt.rc('font', size=14)
plt.rc('font',family='Times New Roman')
plt.rc('font', family='Arial')
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=False, figsize=(8,10))
fig.suptitle('Separator of Li-ion battery')
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
ax1.set_xlim( xmin=0,xmax=3000)
ax2.set_ylim(ymin=0.0, ymax=2500)
ax2.set_xlim(xmin=0, xmax=L)
ax3.set_xlim(xmin=0, xmax=L)
ax1.set_xlabel(r'Salt concentration mol/m3')
ax1.set_ylabel(r'Electrolyte conductivity, S/m')
ax2.set_xlabel(r'Distance microns')
ax2.set_ylabel(r'salt concentration')
ax3.set_ylabel(r'$ \phi_2 $')



x = np.linspace(0, L, num=100)
k_eff = np.zeros(100)
c = np.linspace(0,3000, num=100)
k_eff = conductivity(c)
ax1.plot(c, k_eff, label=r'conductivity')
k_eff = eps**1.5 * conductivity(c)
c = np.zeros(100)



c = ((cd/F)*(1-t_plus)/D_e)*(L/2-x)*1e-6 + c0
k_eff = conductivity(c)

y0 = 0
sol = odeint(dy_dx, y0, x)
ax2.plot(x, c)
ax3.plot(x, sol)
#ax3.plot(x, 1000*sol)