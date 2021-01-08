import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, R, N_A
from scipy import special
import matplotlib.animation as animation
from matplotlib.widgets import Button
#import matplotlib

"""This model (Python 3.8) assumes equilbrium(infinite kinetics at the electrode
A one electron transfer and that diffusion coefficients are equal.  Best
suited to exploring the effect of scan rate and diffusivity have on CV.
The model also includes double layer charging. December 26, 2020
The three parameters you may want to change are sr, the scan rate, C the
capacitance, and D, the diffusion coefficient"""

#TEMPERATURE AND CONSTANTS
T_ex = 298.15   #Temperature, K
F =e*N_A #Faraday's constant
frt = F/R/T_ex

D = 1.0e-9      #Diffusion coefficient, m2/s
L = 0.001 #distance over which plotted
N = 500 #number of time steps

#TIME INTEGRATION PARAMETERS
t_start = 0.01                 #always start at time zero
t_final = 20 # follow for 20 seconds
#t = np.linspace(t_start, t_final, num=N)
t = np.logspace(-1.1,2,N) #1.3 corresponds to 20 seconds

#INITIAL CONCENTRATION
conc = np.zeros((1000, len(t)))
c0= 100 #intial concentration of 100 mol/m^3

x = np.linspace(0,L,1000) #distance from electrode in meters

currentdensity_init = np.zeros(len(t))
currentdensity = np.zeros(len(t))
currentdensity[1:] = F*np.sqrt(D)*c0/np.sqrt(t[1:]*np.pi)



for ii in range (0,1000):
    conc[ii,0] = c0
    
for i in range (0, 1000):
    for j in range(1,N):
        eta = x[i]/np.sqrt(4*D*t[j])
        conc[i, j] = c0*special.erf(eta)

plt.rc('font', size=14)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,6))
line1, = ax1.plot(1000*x, conc[:,[0]], 'k--')
line2, = ax2.plot(t, currentdensity_init, 'b', linewidth='2')

#####################
#adding inset for detailed look at gradient at surface
y_ins = np.zeros(30)
x_ins = np.zeros(30)

y_ins = conc[0:30,0]
x_ins = 1000*x[0:30]

left, bottom, width, height = [0.27, 0.2, 0.2, 0.25]
ax3 = fig.add_axes([left, bottom, width, height])
ax3.set_yticklabels([0,50,100],fontsize=10)
ax3.set_xticklabels([0,0.01,0.02], fontsize=10)
ax3.text(0.005, 87, r'Inset highlighting', fontsize=10)
ax3.text(0.005, 75, r'gradient at surface', fontsize=10)

line3, = ax3.plot(x_ins,y_ins,'g--')

#######################################



stepchange = plt.axes([0.75, 0.80, 0.10, 0.05])
initiate = Button(stepchange, 'initiate', color='dodgerblue', hovercolor='0.975')

pausebutton = plt.axes([0.73, 0.72, 0.15, 0.05])
pause_resuem = Button(pausebutton, 'pause/resume', color='thistle', hovercolor='0.975')

#######################   ANIMATION
def run_animation():
    anim_running = True

    def onClick(event):
        nonlocal anim_running
        if anim_running:
            anim.event_source.stop()
            anim_running = False
        else:
            anim.event_source.start()
            anim_running = True

    def animFunc(k):
        line1.set_ydata(conc[:,[k]])
        line2.set_ydata(currentdensity[1:k])
        line2.set_xdata(t[1:k])
        line3.set_ydata(conc[0:30, [k]])# update the data.
        # Animation update function here

    fig.canvas.mpl_connect('button_press_event', onClick)

    anim = animation.FuncAnimation(fig, animFunc, frames=499, repeat=False, interval=3)

#run_animation()


def potentialstep(event):
    run_animation()

initiate.on_clicked(potentialstep)  



fig.suptitle('Potential step experiment', fontsize=18, y=0.95)
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim([0,120])
ax1.set_xlim([0,1000*L])
ax2.set_xlim([-1, 100])
ax2.set_ylim([0,500])
ax1.set_ylabel(r'Concentration, $\mathrm {mol\/m^{-3}}$')
ax1.set_xlabel(r'Distance from electrode, mm')
ax2.set_xlabel(r'Time, s')
ax2.yaxis.set_label_position("right")
ax2.set_ylabel(r'Current density, $\mathrm {A\/m^{-2}}$')
ax3.set_ylim([0,120])
ax3.set_xlim([0, 0.025])

anno1 = ax1.annotate('', xytext=(0,0), xy=(0.42,14), arrowprops=dict(color='gainsboro', arrowstyle="-"))
anno2 = ax1.annotate('', xytext=(0,0.025), xy=(0.99,14), arrowprops=dict(color='gainsboro', arrowstyle="-"))
plt.minorticks_on()

plt.show()
################################






























        



