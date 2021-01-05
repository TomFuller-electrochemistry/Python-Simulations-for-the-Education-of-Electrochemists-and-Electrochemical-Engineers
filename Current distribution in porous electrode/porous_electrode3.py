"""
In this example, a slider is used to choose the ratio of the ionic and electronic
conductivity as well as the Wa for a porous electrode with linear kinetics."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

fig, ax = plt.subplots(figsize=(8, 8))
plt.subplots_adjust(left=0.15, bottom=0.25)
ax.text(0.25, 3.5, r'$\nu^2$ is defined by Eq. 5.31', fontsize=14)
ax.text(0.25, 3.2, 'similar to inverse of Wa', fontsize=14)
ax.text(0.02, 4.2, 'Current distribution across a porous electrode', fontsize=18)
ax.text(-0.1, 1.0, 'Dimensionless reaction rate', fontsize=16, rotation=90)
ax.text(0.1, 5.1, 'Dimensionless distance from current collector', fontsize=16)
plt.xlim(0.0,1.0)
plt.ylim(0.0,5.0)
plt.grid(True)
t = np.arange(0.0, 1.0, 0.0001)#dimensionless distance from current collector
Kr0 = 1.000
nu0 = 3.000
delta_x = 0.001

s = (nu0*np.cosh(nu0*t) + nu0*Kr0*np.cosh(nu0*(t-1.000)))/((1.000+Kr0)*np.sinh(nu0))
l, = plt.plot(t, s, lw=2, color='k')

axnu = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor='aliceblue')
axKr = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor='linen')
ax.text(0.25, -1.35, 'ratio of ionic and electronic conductivity', fontsize=10)
ax.text(0.30, -0.5, 'ohmic vs kinetic resistance', fontsize=10)

sKr = Slider(axKr, r'$K_r = \kappa/ \sigma$', 0.01, 30.0, valinit=Kr0, color='peru', valstep=delta_x)
snu = Slider(axnu, r'$\nu$', 0.01, 10.0, valinit=nu0, color='dodgerblue')

def update(val):
    nu = snu.val
    Kr = sKr.val
    l.set_ydata((nu*np.cosh(nu*t) + nu*Kr*np.cosh(nu*(t-1.000)))/((1.000+Kr)*np.sinh(nu)))
    fig.canvas.draw_idle()

sKr.on_changed(update)
snu.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color='silver', hovercolor='0.975')

def reset(event):
    sKr.reset()
    snu.reset()
button.on_clicked(reset)
    
plt.show()
