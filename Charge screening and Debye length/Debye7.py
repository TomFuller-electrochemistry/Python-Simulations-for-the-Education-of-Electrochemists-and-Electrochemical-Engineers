import numpy as np
from scipy.constants import epsilon_0, e, R, N_A
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.ticker import (MultipleLocator)

#SET UP PLOT
plt.rc('font', size=14)
plt.rc('font',family='Times New Roman')
plt.rc('font', family='Arial')
fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(8,8))
fig.suptitle('Charge Screening', fontsize=18, y=0.93)
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim( ymin=0, ymax=0.2)
ax2.set_ylim(ymin=0.0, ymax=2.0)
plt.xlim([0,10])
ax1.xaxis.set_major_locator(MultipleLocator(1))
ax2.set_xlabel(r'Distance from central positive ion, nm')
ax1.set_ylabel(r'Potential, V')
ax2.set_ylabel(r'$ c/c^{\infty}$')
ax1.text(7.05, 0.16, r'$\phi = \frac{z_c q}{4\pi\epsilon_0 r}$', fontsize=14)

#ESTABLISH TEMPERATURE AND PARAMETERS
T = 298 #temperature, K
F = e*N_A #Faraday's constant
ce = 200 #default salt concentration. mol/m^3
a = 0.2e-9 #assumed size of charge, nm
rp = 78.3 #relative permittivity of water at 25 C
qT = 1 #charge on ion

#SET UP SLIDER BARS TO ADJUST THE CONCENTRATION OF ELECTROLYTE
#corse control
ax2nu = plt.axes([0.38, 0.20, 0.35, 0.02], facecolor='aliceblue')
ion_strength = Slider(ax2nu, r'$c_e$', 2, 990, valinit=ce, color='dodgerblue')
ion_strength.valtext.set_visible(False)

#fine control
ax2f_c = plt.axes([0.48, 0.14, 0.15, 0.01], facecolor='aliceblue')
fine_control = Slider(ax2f_c, r'fine', 0, 10, valinit=0, color='dodgerblue')
fine_control.valtext.set_visible(False)

ax2.text(3.7, 0.35, r'$\mathrm {electrolyte\/concentration, mol/m^3}$', fontsize=12)

#add text box that shows the electrolyte concentration
electrolyte = ax2.text(8.3, 0.35, '',bbox=dict(facecolor='aliceblue'))
electrolyte.set_text(r'$c_e$ = %.1f'%  ce)


#DEFINE FUNCTION USED
def calc_debye_length(ce):
    lamde = np.sqrt(epsilon_0*rp*R*T/F/F/ce)
    return lamde

def calc_unscreened_potential(r, qT):
    return qT * e / 4 / np.pi / (epsilon_0*rp) / r

def calc_e_potential(r, lam_De, qT):
    return calc_unscreened_potential(r, qT) * np.exp((a-r) / lam_De)/(1+a/lam_De)

def calc_unscreened_vacuum(r, qT):
    return qT * e / 4 / np.pi / epsilon_0 / r

def update(val):
    course_c = ion_strength.val
    cf = fine_control.val
    ce = course_c + cf
    electrolyte.set_text(r'$c_e$ = %.1f'%  ce)
    lam_De = calc_debye_length(ce)
    phi = calc_e_potential(r, lam_De, qT)
    c_an = np.exp(F*phi/R/T)
    c_cat = np.exp(-F*phi/R/T)
    ca.set_ydata(c_an)
    cc.set_ydata(c_cat)
    ps.set_ydata(phi)
    xd = lam_De*1e9
    ld1.set_xdata(xd)
    ld2.set_xdata(xd)
    fig.canvas.draw_idle()
    return



lam_De = calc_debye_length(ce)

# range of distances to plot phi for, in m.
rmin = a
rmax = 1.0e-8
r = np.linspace(rmin, rmax, 100)


phi_unscreened = calc_unscreened_potential(r, qT)
phi = calc_e_potential(r, lam_De, qT)
phi_vac = calc_unscreened_vacuum(r, qT)


c_cat = np.exp(-F*phi/R/T)
c_an = np.exp(F*phi/R/T)

# Plot the figure.
ax1.plot(r*1e9, phi_vac, label=r'vacuum')
ax1.plot(r*1.e9, phi_unscreened, label=r"pure water")
ps, = ax1.plot(r*1.e9, phi, label=r'aqueous electrolyte')

ax1.legend(bbox_to_anchor=(0.65, 0.6))

cc, = ax2.plot(r*1e9, c_cat, label=r'cation')
ca, = ax2.plot(r*1e9, c_an, label=r'anion')
yd = np.linspace(0, 2, num=20)
xd = np.zeros(20)
xd[0:20] = lam_De*1e9
ld1, = ax1.plot(xd, yd, ls='--', c='k')
ld2, = ax2.plot(xd, yd, ls='--', c='k', label=r'Debye length')
ax2.legend(bbox_to_anchor=(0.65, 0.6))

ion_strength.on_changed(update)
fine_control.on_changed(update)
plt.show()