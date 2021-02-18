"""
In this example, sliders are available to explore three parameters used
in a full cell: the conductivity of solution, the exchange-current 
density, and the current density.
The nominal values reproduce Figure 3.10 from the text.
"""
#IMPORT THE NEEDED MODULES
import numpy as np
from scipy.constants import e, R, N_A
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.optimize import fsolve
from matplotlib.patches import Rectangle

#get rid of this later
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

#SET UP THE PLOTTING AREA
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(12,8))
fig.text(0.20, 0.14, "Exchange current density for Zn", fontsize=12, ha='center')
fig.text(0.20, 0.08, "Current density", fontsize=12)
fig.text(0.20, 0.02, "Electrolyte conductivity", fontsize=12)
plt.subplots_adjust(left=0.15, bottom=0.25)
fig.suptitle('Ni-Zn Full cell with Different Reference Electrodes', fontsize=18, y=0.93)
ax1.set_xlabel(r'Distance from Zn Electrode, mm', fontsize=14)
ax1.set_ylabel(r'Potential, V', fontsize=16)
ax1.set_ylim(-0.2,1.95)
ax1.set_xlim(0.0, 2.0)
ax2.set_xlim(0.0, 2.0)
ax1.grid(True, which='both')
ax2.grid(True, which='both')
plt.rc('font', size=14)


#SPECIFY PARAMETERS
T = 298.15 #temperature K
F = e*N_A # Faraday's constant
frt = F/R/T
#negative electrode data
R_aZn =2.0
xi0n = R_aZn*60 #effective exchange current density [A cm^-2]
alpha_an = 1.5
alpha_cn = 0.5
#positive electrode data
R_aNi =100
xi0p = R_aNi*0.61 #effective exchange current density [A cm^-2]
alpha_ap = 0.5
alpha_cp = 0.5

kappa = 60 #conductivity S/m
cd = 1000 #current density, A m-2

#Define functions that use Butler-Volmer kinetics to express current density
def ne(z): return cd - xi0n*(np.exp(alpha_an*frt*z)-np.exp(-alpha_cn*frt*z))
def pe(z2): return cd + xi0p*(np.exp(alpha_ap*frt*z2)-np.exp(-alpha_cp*frt*z2))
#use fsolve to find overpotential    
eta_sn = fsolve(ne, 1)
eta_sp = fsolve(pe, 1)

#arbitrarily set potential of negative electrode to zero

#calculate solution potential assuming a Zn reference electrode
phi_2 = np.zeros(50)
x = np.linspace(0,2,num=50) #assumes thickness of 2 mm
phi_2 = -eta_sn - cd*x/kappa/1000

#ohmic drop
eta_Omega = phi_2[0] - phi_2[49]
#convert to HgO reference (Un=0.861)
Un = 0.861
phi_2HgO = np.zeros(50)
phi_2HgO = phi_2 + Un

#calculate cell potential
phi_1p = 1.74 + eta_sp + phi_2[49]

line1, = ax1.plot(x, phi_2, label=r'$\phi_2$ (Zn ref)')
y_eta_sn = np.linspace(0, -eta_sn, num=10)
xd = np.zeros(10)
xd[0:10] = 0.03
line2, = ax1.plot(xd, y_eta_sn, ls='--', c='k')

line3, = ax2.plot(x, phi_2HgO, label=r'$\phi_2$ (HgO ref)')

#y_eta_sp = np.linspace(1.74+phi_2[49], 1.74+eta_sp+phi_2[49], num=10)
#xp = np.zeros(10)
#xp[0:10] = 1.97
#line3, = ax.plot(xp, y_eta_sp, ls='--', c='k')

#add arrows to point out overpotentials and metal potentials
anno1 = ax1.annotate('$ \eta_{s,n}$', xytext=(0.25,0.125), xy=(0.03,-eta_sn/2), arrowprops=dict(arrowstyle="->"), fontsize=14)
anno2 = ax1.annotate('$ \eta_{s,p}$', xytext=(1.6,1.5), xy=(1.26,1.74-eta_sn+eta_sp/2), arrowprops=dict(arrowstyle="->"), fontsize=14)
anno3 = ax1.annotate(r'$ \phi_{1,n} =0$', xytext=(0.1, 0.3,), xy=(0,0), arrowprops=dict(arrowstyle="->"), fontsize=14)
anno4 = ax1.annotate(r'$ \phi_{1,p}$', xytext=(1.6, 1.75,), xy=(2,phi_1p), arrowprops=dict(arrowstyle="->"), fontsize=14)
anno5 = ax2.annotate(r'$ \phi_{1,p}$', xytext=(1.6, 1.75,), xy=(2,phi_1p), arrowprops=dict(arrowstyle="->"), fontsize=14)
anno6 = ax2.annotate(r'$ \phi_{1,n} =0$', xytext=(0.1, 0.3,), xy=(0,0), arrowprops=dict(arrowstyle="->"), fontsize=14)

ax1.text(0.865, 0.54, r'$ U_{cell} =1.74$', fontsize=14, rotation=90, color='w')
ax2.text(0.1, 1.8, r'HgO/Hg reference', fontsize=14, color='k')
ax1.text(0.1, 1.8, r'Zn reference', fontsize=14, color='k')

blku = ax1.add_patch(Rectangle((0.75, -eta_sn), 0.25, 1.74, facecolor='k'))
blketa_n = ax1.add_patch(Rectangle((1.0, -eta_sn), 0.25, eta_sn, facecolor = 'dodgerblue'))
blketa_p = ax1.add_patch(Rectangle((1.0, 1.74-eta_sn+eta_sp), 0.25, -eta_sp, facecolor = 'blue'))
blkW = ax1.add_patch(Rectangle((1.0, 0), 0.25, phi_2[0]-phi_2[49], facecolor = 'peru'))

cell_pot = ax1.text(1.26, 0.9, '',bbox=dict(fc='aliceblue', lw='1'))
cell_pot.set_text(r'$V _{cell}$ = %.3f'%  phi_1p)
overpotentialn = fig.text(0.8,0.15, '',bbox=dict(fc='aliceblue', lw='1' ))
overpotentialn.set_text(r'$\eta_{s,\mathrm{Zn}}$ = %.3f'%  eta_sn)
overpotentialp = fig.text(0.8,0.09, '',bbox=dict(fc='aliceblue', lw='1' ))
overpotentialp.set_text(r'$\eta_{s,\mathrm{Ni}}$ = %.3f'%  eta_sp)
overpotentialW = fig.text(0.8,0.03, '',bbox=dict(fc='aliceblue', lw='1' ))
overpotentialW.set_text(r'$\eta_{\Omega}}$ = %.3f'%  eta_Omega)

#SET UP SLIDER BARS
axXchange = plt.axes([0.4, 0.14, 0.30, 0.02], facecolor='aliceblue')
axcd = plt.axes([0.4, 0.08, 0.30, 0.02], facecolor='linen')
axkappa = plt.axes([0.4, 0.02, 0.30, 0.02], facecolor='violet')
skappa = Slider(axkappa, r'$ \kappa$', 30, 200, valinit=60, color='peru', valstep=0.001, valfmt='%1.1f')
sXchange = Slider(axXchange, r'$i_{o\mathrm{, Zn}}$', 1, 300, valinit=60, color='dodgerblue')
scd = Slider(axcd, r'$i$', -500, 1500, valinit=1000, color='indigo')

#Function to update plot when sliders are moved
def update(val):
    xi0n = R_aZn*sXchange.val
    cd = scd.val
    kappa = skappa.val
    def ne(z2): return cd - xi0n*(np.exp(alpha_an*frt*z2)-np.exp(-alpha_cn*frt*z2))
    eta_sn = fsolve(ne, 1)
    def pe(z3): return cd + xi0p*(np.exp(alpha_ap*frt*z3)-np.exp(-alpha_cp*frt*z3))
    eta_sp = fsolve(pe,1)
    phi_2 = -eta_sn - cd*x/kappa/1000 #1000 is conversion factor
    phi_2HgO = phi_2 + Un
    line1.set_ydata(phi_2)
    line3.set_ydata(phi_2HgO)
    new_overpotentialn = np.linspace(0, -eta_sn, num=10)
    line2.set_ydata(new_overpotentialn)
    anno1.xy = (0.03,-eta_sn/2)
    anno2.xy = (1.26, 1.74-eta_sn+eta_sp/2)
    blku.xy = (0.75,-eta_sn)
    blketa_n.xy = (1.0, -eta_sn)
    blketa_n.set_height(eta_sn)
    blketa_p.xy = (1.0, 1.74-eta_sn+eta_sp)
    blketa_p.set_height(-eta_sp)
    blkW.set_height(phi_2[0]-phi_2[49])
    phi_1p = 1.74 + eta_sp + phi_2[49]
    anno4.xy = (2, phi_1p)
    anno5.xy = (2, phi_1p)
    eta_Omega = abs(phi_2[0] - phi_2[49])
    cell_pot.set_text(r'$V _{cell}$ = %.3f'%  phi_1p)
    overpotentialn.set_text(r'$\eta_{s,\mathrm{Zn}}$ = %.3f'%  eta_sn)
    overpotentialp.set_text(r'$\eta_{s,\mathrm{Ni}}$ = %.3f'%  eta_sp)
    overpotentialW.set_text(r'$\eta_{\Omega}}$ = %.3f'%  eta_Omega)
    #overpotential.set_text(r' = %.3f'% eta_sn)
    
    fig.canvas.draw_idle()

sXchange.on_changed(update)
scd.on_changed(update)
skappa.on_changed(update)

ax1.legend(bbox_to_anchor=(0.55, 0.25))  
ax2.legend(bbox_to_anchor=(0.7, 0.4))    
plt.show()
