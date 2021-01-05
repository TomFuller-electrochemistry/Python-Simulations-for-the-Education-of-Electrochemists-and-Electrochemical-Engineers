# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 20:10:20 2020

@author: tomff
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
import matplotlib.image as mpimg#from matplotlib import animation as an
import control
from matplotlib.widgets import Slider

# Parameters for the circuit
sigma = 5
R_ohm = 1 #solution resistance ohm-cm2
R_f = 5 #charge transfer resistance ohm-cm2
C = 0.00005 #double layer capacitance F/cm2

img = mpimg.imread('Randlescircuit.png')#get image for plot

# Set up the figure
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8,8))
plt.subplots_adjust(hspace=0.2)
plt.subplots_adjust(wspace=0.3)
ax3.sharex(ax1)
ax4.sharex(ax2)

ax1.grid(True)
ax2.grid(True, which='both')
ax4.grid(True)
ax1.set_xlim(-0.1, 15)
ax1.set_ylim(0, 15)
ax2.set_xlim(0.01, 100000)
ax2.set_ylim(0.1,100)
ax2.set_yscale("log")
ax3.xaxis.set_visible(False)
ax3.yaxis.set_visible(False)
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)
ax3.spines["left"].set_visible(False)
ax3.spines["bottom"].set_visible(False)
ax3.text(2.0, 0.1, r'Adjustable Parameters', fontsize=12)
ax2.text(1, 1.1, r'Bode plot', fontsize=14)
ax1.text(1, 6.5, r'Nyquist plot', fontsize=14)
ax4.set_ylim(-90, 0)
ax4.set_xscale("log")


w = np.logspace(-1.2,6.8,1000)

# calculate impedance, and make Nyquist plot
a1 = 1 + sigma*np.sqrt(w)*C
b1 = w*C*(sigma/np.sqrt(w)+R_f)
c1 = sigma/np.sqrt(w) + R_f
d1 = -sigma/np.sqrt(w)

Wr = R_ohm + (a1*c1+b1*d1)/(a1*a1+b1*b1)
Wi = -(a1*d1-c1*b1)/(a1*a1+b1*b1)
#Wr = R_ohm + R_f/(1+w*w*R_f*R_f*C*C) + sigma/np.sqrt(w)
#Wi = w*R_f*R_f*C/(1+w*w*R_f*R_f*C*C)  + sigma/np.sqrt(w)
ax1.set_ylabel(r'-Zi, $\Omega \mathrm {-cm^2}$')
ax1.set_xlabel(r'Zr, $\Omega \mathrm {-cm^2}$')
ax2.set_ylabel('Amplitude ratio')
ax4.set_xlabel('frequency, Hz')
ax4.set_ylabel('Phase angle, degrees')
p1, = ax1.plot(Wr, Wi, color='k')
anno1 = ax1.annotate('1 kHz', xytext=(1,4), xy=(Wr[625],Wi[625]), arrowprops=dict(arrowstyle="->"))
anno2 = ax1.annotate('1 Hz', xytext=(5,5), xy=(Wr[150],Wi[150]), arrowprops=dict(arrowstyle="->"))

# calculate gain and phase angle, make Bode plot
AR = np.sqrt(Wr*Wr + Wi*Wi)
p2, = ax2.plot(0.159*w, AR)

PAng = 57.3*np.arctan(-Wi/Wr)
p4, = ax4.plot(0.159*w, PAng)

#SET UP SLIDER BARS
ax3rw = plt.axes([0.11, 0.2, 0.25, 0.02], facecolor='aliceblue')
ohmres = Slider(ax3rw, r'$R_{\Omega}$', 0.01, 10, valinit=1, color='dodgerblue', valfmt='%1.2f')
ax1up = plt.axes([0.11, 0.25, 0.25, 0.02], facecolor='aliceblue')
capacitance = Slider(ax1up, r'$C$', 1e-5, 1e-4, valinit=5e-5, color='dodgerblue', valfmt='%1.6f')
ax1Rf = plt.axes([0.11, 0.30, 0.25, 0.02], facecolor='aliceblue')
kineticres = Slider(ax1Rf, r'$R_f$', 0.001, 10, valinit=5, color='dodgerblue', valfmt='%1.2f')
ax1Warb = plt.axes([0.11, 0.35, 0.25, 0.02], facecolor='aliceblue')
Warburg = Slider(ax1Warb, r'Warburg', 0.01, 10, valinit=5, color='dodgerblue')

def update(val):
    R_ohm = ohmres.val
    sigma = Warburg.val
    R_f = kineticres.val
    C = capacitance.val
    a1 = 1 + sigma*np.sqrt(w)*C
    b1 = w*C*(sigma/np.sqrt(w)+R_f)
    c1 = sigma/np.sqrt(w) + R_f
    d1 = -sigma/np.sqrt(w)

    Wr = R_ohm + (a1*c1+b1*d1)/(a1*a1+b1*b1)
    Wi = -(a1*d1-c1*b1)/(a1*a1+b1*b1)
    #Wr = R_ohm + R_f/(1+w*w*R_f*R_f*C*C) + sigma/np.sqrt(w)
    #Wi = w*R_f*R_f*C/(1+w*w*R_f*R_f*C*C)  + sigma/np.sqrt(w)
    p1.set_ydata(Wi)
    p1.set_xdata(Wr)
    AR = np.sqrt(Wr*Wr + Wi*Wi)
    p2.set_ydata(AR)
    PAng = 57.3*np.arctan(-Wi/Wr)
    p4.set_ydata(PAng)
#   anno1.set_position((Wr[625],Wi[625]))
    anno1.xy = (Wr[625],Wi[625])  
    anno2.xy = (Wr[150],Wi[150])  
#   anno1 = ax1.annotate('1 kHz', xytext=(1,4), xy=(Wr[625],Wi[625]), arrowprops=dict(arrowstyle="->"))
#   anno2 = ax1.annotate('1 Hz', xytext=(5,5), xy=(Wr[150],Wi[150]), arrowprops=dict(arrowstyle="->"))

    fig.canvas.draw_idle()
    
    
# add image for Randles Circuit
imagebox_python = OffsetImage(img, zoom=.45)
xy = [5.5,11.7]
ab_RandlesLogo = AnnotationBbox(imagebox_python, xy, xybox=(1., -1.), boxcoords='offset points')
ax1.add_artist(ab_RandlesLogo)

fig.suptitle('EIS for Randles Circuit', fontsize=18, y=0.93)
ohmres.on_changed(update)
capacitance.on_changed(update)
kineticres.on_changed(update)
Warburg.on_changed(update)
plt.show()
