"""
In this example, sliders are available to explore three parameters used in the
creation of a polization curve: the ohmic resistance, the exchange-current 
density, and the limiting current.
The nominal values are those of problem 9.4 from the text.
"""
#IMPORT THE NEEDED MODULES
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.optimize import fsolve

#SET UP THE PLOTTING AREA
fig, ax = plt.subplots(figsize=(8, 8))
plt.subplots_adjust(left=0.15, bottom=0.25)
ax.text(0.75, 0.23, 'Current density, A-$\mathrm{cm^{-2}}$', fontsize=16)
ax.text(-0.30, 0.50, 'Cell potential, V', fontsize=18, rotation=90)
plt.xlim(0.0,3.0)
plt.ylim(0.3,1.0)
plt.grid(True)

#SPECIFY PARAMETERS
delta_x = 0.001
xi0 = 0.020 #exchange current density [A cm^-2]
R0 = 0.10 #resistance ohm-cm^2
ilim0 = 1.4 #limiting current density, A m^-2
nk = int(3.0/delta_x)
cd = np.arange(0.0001, 3.0, delta_x)#current density, A m^-2

###################DEFINE FUNCTIONS
def ys(n1, k1, mt, cdm): #polarization
    def fk(xk, *data): #Butler-Volmer equation
        p1, p2 = data
        frt = 19.462 # alpha*F/RT
        return p2 - p1*(np.exp(frt*xk)-np.exp(-frt*xk))
    
    xk0 = 0.5
    eta_s = np.arange(0, 3.0, 0.001)
    mtp   = np.arange(0,3.0, 0.001)
    for ik1 in range(0,nk):
            data = (k1, cdm[ik1])
            eta_t = fsolve(fk, xk0, args=data) #use built-in solver
            eta_s[ik1] = eta_t
            if cdm[ik1] < mt:
                mtp[ik1] = cdm[ik1]/(mt - cdm[ik1])
            else:
                mtp[ik1] = 50
    ts = 0.0591 # Tafel slope of 59 mV/decade for mass-transfer
    ys = 1.05 - n1*cdm - eta_s - ts*mtp
#   ys = 1.05 - n1*cdm - eta_s - ts*(cdm/(mt-cdm))
    return (ys)
####################    

#PLOT INITIAL LINE
w0 = ys(R0, xi0, ilim0, cd)# resistance, current density
mp, = plt.plot(cd, w0, lw=2, color='k')

#SET UP SLIDER BARS
axR = plt.axes([0.35, 0.14, 0.50, 0.03], facecolor='aliceblue')
axxi = plt.axes([0.35, 0.08, 0.50, 0.03], facecolor='linen')
axmt = plt.axes([0.35, 0.83, 0.50, 0.03], facecolor='violet')

ax.text(1.2, 0.08, 'exchange current density, A-$\mathrm{cm^{-2}}$', fontsize=10)
ax.text(1.2, 0.15, 'ohmic resistance, $\Omega$-$\mathrm{cm^{2}}$', fontsize=10)
ax.text(1.2, 0.91, 'limiting current density, A-$\mathrm{cm^{-2}}$', fontsize=10)

sxi = Slider(axxi, r'0.0001', 0.0001, 0.5, valinit=xi0, color='peru', valstep=0.001, valfmt='%1.4f')
sR = Slider(axR, r'0.01', 0.003, 3.0, valinit=R0, color='dodgerblue')
smt = Slider(axmt, r'1.0', 1.01, 3.0, valinit=ilim0, color='indigo')

#UPDATE MODEL
def update(val):
    R = sR.val
    Kr = sxi.val
    mt = smt.val
    mp.set_ydata(ys(R, Kr, mt, cd))
    fig.canvas.draw_idle()

sxi.on_changed(update)
sR.on_changed(update)
smt.on_changed(update)

#RESET TO INITIAL VALUES
resetax = plt.axes([0.75, 0.70, 0.1, 0.04])
button = Button(resetax, 'Reset', color='silver', hovercolor='0.975')

def reset(event):
    sxi.reset()
    sR.reset()
    smt.reset()
    
button.on_clicked(reset)

    
plt.show()
