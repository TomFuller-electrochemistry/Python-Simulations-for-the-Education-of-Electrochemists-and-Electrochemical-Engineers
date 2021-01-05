# -*- coding: utf-8 -*-
"""
Created on Mon Dec  26 10:40:57 2020

@author: tomff
"""
import numpy as np
#from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.constants import e, R, N_A
from matplotlib.widgets import Slider

L = 20 #thickness of separator, microns
cd = 100 #current density, A/m^2
F = e*N_A #Faraday's constant, C/eq
#R = 8.314 #universal gas constant
T = 298 # temperature, K
c0 = 250 #concentration of electrolyte, mol/m^3
eps = 0.3
#t_plus = 0.19 #lithium ion transference number
#D_b = 2.79e-10 #diffusivity, m^2/s
#D_e = eps**1.5 * D_b
#mobility data
u_p = 2.8e-13
u_n = 4.3e-13

D = 2*R*T*u_p*u_n/(u_p+u_n) #m^2/s
tp = u_p/(u_p+u_n)
kappa = F*F*c0*(u_p+u_n) # S/m

#SET UP PLOT
plt.rc('font', size=14)
plt.rc('font',family='Times New Roman')
plt.rc('font', family='Arial')
fig, (ax1, ax2) = plt.subplots(2, sharex=False, figsize=(8,8))
fig.suptitle('Dilute binary 1:1 electrolyte', fontsize=18, y=0.93)
ax1.grid(True)
ax2.grid(True)
ax1.set_xlim( xmin=0,xmax=500)
ax1.set_ylim(ymin=0, ymax=5.0)
ax2.set_ylim(ymin=200, ymax=300)
ax2.set_xlim(xmin=0, xmax=L)
ax1.set_xlabel(r'Salt concentration, mol $\mathrm {m^{-3}}$')
ax1.set_ylabel(r'Electrolyte conductivity, S/m')
ax2.set_xlabel(r'Distance, $\mathrm {\mu  m}$')
ax2.set_ylabel(r'salt concentration, mol $\mathrm {m^{-3}}$')
#ax3.set_ylabel(r'$ \phi_2 $')

#add cell potential text box
trans_num = ax1.text(300, 1.1, '',bbox=dict(facecolor='aliceblue'))
trans_num.set_text(r'$t _{+}$ = %.2f'%  tp)
diff_coef = ax1.text(300, 0.5, '',bbox=dict(facecolor='aliceblue'))
diff_coef.set_text(r'$D \mathrm {[ m^2 s^{-1}]}$ = %.1e'%  D)

#SET UP SLIDER BARS
ax1up = plt.axes([0.165, 0.8, 0.25, 0.02], facecolor='aliceblue')
mobilityp = Slider(ax1up, r'$u_{+}$', 0.01, 10, valinit=5, color='dodgerblue', valfmt='%1.1f $\mathrm { X10^{13}m^2 mol\/ J^{-1} s^{-1}}$')
ax1un = plt.axes([0.165, 0.75, 0.25, 0.02], facecolor='aliceblue')
mobilityn = Slider(ax1un, r'$u_{-}$', 0.01, 10, valinit=5, color='dodgerblue', valfmt='%1.1f $\mathrm { X10^{13}m^2 mol\/ J^{-1} s^{-1}}$')
#ax2.text(3.7, 0.35, r'$\mathrm {electrolyte\/concentration, mol/m^3}$', fontsize=12)

#ax2c = plt.axes([0.20, 0.2, 0.25, 0.02], facecolor='aliceblue')
#concentration = Slider(ax2c, r'$c_e$', 1, 250, valinit=c0, color='dodgerblue')
ax2cd = plt.axes([0.40, 0.40, 0.25, 0.02], facecolor='aliceblue')
currentdensity = Slider(ax2cd, r'$i$', -500, 500, valinit=100, color='dodgerblue')
ax2.text(8, 275, r'$\mathrm {current\/ density, A/m^2}$', fontsize=12)

def conductivity(c, u_p, u_n):
    k_eff = F*F*c*(u_p+u_n) # S/m
    return k_eff

def update(val):
    mp = mobilityp.val
    mn = mobilityn.val
    cd = currentdensity.val
    u_p = mp*1e-13
    u_n = mn*1e-13
    c_an = conductivity(c,u_p, u_n)
    k1.set_ydata(c_an)
    tp = u_p/(u_p+u_n)
    trans_num.set_text(r'$t_{+}$ = %.2f'%  tp)
    D = 2*R*T*u_p*u_n/(u_p+u_n) #m^2/s
    diff_coef.set_text(r'$D[m^2 s^{-1}]$ = %.2e'%  D)
    c_cur = c0 - (cd/2/F/R/T/u_p)*(x-L/2)*1e-6
    cdp.set_ydata(c_cur)
    fig.canvas.draw_idle()
    
def update2(val):
    cd = currentdensity.val
    mp = mobilityp.val
    u_p = mp*1e-13
    c_cur = c0 - (cd/2/F/R/T/u_p)*(x-L/2)*1e-6
    cdp.set_ydata(c_cur)
    fig.canvas.draw_idle()

x = np.linspace(0, L, num=100)
k_eff = np.zeros(100)
c = np.linspace(0,500, num=100)
k_eff = conductivity(c, u_p, u_n)
k1, = ax1.plot(c, k_eff, label=r'conductivity')
#k_eff = eps**1.5 * conductivity(c, u_p, u_n)


c_cur = c0 - (cd/2/F/R/T/u_p)*(x-L/2)*1e-6
k_eff = conductivity(c, u_p, u_n)

y0 = 0
#sol = odeint(dy_dx, y0, x)
cdp, = ax2.plot(x, c_cur)
#ax3.plot(x, sol)
#ax3.plot(x, 1000*sol)
mobilityp.on_changed(update)
mobilityn.on_changed(update)
currentdensity.on_changed(update)
plt.show()