# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 20:10:20 2020

@author: tomff
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import e, R, N_A
#from matplotlib.offsetbox import AnnotationBbox, OffsetImage
#import matplotlib.image as mpimg#from matplotlib import animation as an
from matplotlib.widgets import Slider

# Parameters for the circuit
xcur = 0.1 #A/m^2
F = e*N_A #Faraday's constant
T = 298 #Temperature, K
#R = 8.314 #Gas constant, J/mol-K
alpha_a = 0.5
alpha_c = 1.0 - alpha_a

xf = alpha_a*F/R/T

eta_s = np.linspace(0,0.2,1000)

#img = mpimg.imread('Randlescircuit.png')#get image for plot

# Set up the figure
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,6))
ax1.grid(True, which='both')
ax1.set_xlim(1e-3, 10)
ax1.set_ylim(0, 0.2)
ax1.set_xscale("log")
ax1.text(1, 6.5, r'Nyquist plot', fontsize=14)
ax1.set_ylabel(r'Surface overpotential, V', fontsize=14)
ax1.set_xlabel(r'Current density, $\mathrm {A-m^{-2}}$', fontsize=14)
ax2.set_xlim(0,1)
ax2.set_xlabel(r'Current density, $\mathrm {A-m^{-2}}$', fontsize=14)
ax2.grid(True, which='both')

# calculate current density
cd_Taf = xcur*np.exp(xf*eta_s)
p1, = ax1.plot(cd_Taf, eta_s, color='r') #Tafel plot
p1b, = ax2.plot(cd_Taf, eta_s, color='r', label=r'Tafel kinetics')

cd_lin = xcur*(alpha_a+alpha_c)*F*eta_s/R/T
p2, = ax1.plot(cd_lin, eta_s, color='b') #Linear approximation
p2b, = ax2.plot(cd_lin, eta_s, color='b', label=r'Linear kinetics')

cd_bv = xcur*(np.exp(xf*eta_s)-np.exp(-xf*eta_s))
p3, = ax1.plot(cd_bv, eta_s, color='k')
p3b, = ax2.plot(cd_bv, eta_s, color='k', label=r'Butler_Volmer')

ax2.legend(bbox_to_anchor=(0.5, 0.6))
#SET UP SLIDER BARS
ax2xcur = plt.axes([0.65, 0.2, 0.20, 0.02], facecolor='maroon')
exchangecurrentdensity = Slider(ax2xcur, r'log $i_o$', -3, 1, valinit=-1, color='teal', valfmt='%1.2f')
ax1tc = plt.axes([0.15, 0.65, 0.15, 0.02], facecolor='maroon')
trcoef = Slider(ax1tc, r'$\alpha_a$', 0, 1, valinit=0.5, color='teal', valfmt='%1.2f')
#ax1Rf = plt.axes([0.11, 0.30, 0.25, 0.02], facecolor='aliceblue')
#kineticres = Slider(ax1Rf, r'$R_f$', 0.001, 10, valinit=5, color='dodgerblue', valfmt='%1.2f')


def update(val):
    alpha_a=trcoef.val
    alpha_c=1-alpha_a
    xf = alpha_a*F/R/T
    logio = exchangecurrentdensity.val
    xcur = 10**logio
    cd_Taf = xcur*np.exp(xf*eta_s)
    p1.set_xdata(cd_Taf)
    p1b.set_xdata(cd_Taf)
    cd_lin = xcur*(alpha_a+alpha_c)*F*eta_s/R/T
    p2.set_xdata(cd_lin)
    p2b.set_xdata(cd_lin)
    cd_bv = xcur*(np.exp(xf*eta_s)-np.exp(-xf*eta_s))
    p3.set_xdata(cd_bv)
    p3b.set_xdata(cd_bv)
    fig.canvas.draw_idle()
    
    
# add image for Randles Circuit
#imagebox_python = OffsetImage(img, zoom=.45)
#xy = [5.5,11.7]
#ab_RandlesLogo = AnnotationBbox(imagebox_python, xy, xybox=(1., -1.), boxcoords='offset points')
#ax1.add_artist(ab_RandlesLogo)

fig.suptitle('Simplified Forms of the BV Equation', fontsize=18, y=0.93)
exchangecurrentdensity.on_changed(update)
trcoef.on_changed(update)
plt.show()
