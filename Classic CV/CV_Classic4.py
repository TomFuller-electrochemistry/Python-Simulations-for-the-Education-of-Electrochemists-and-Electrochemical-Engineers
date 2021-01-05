import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, R, N_A
import matplotlib.animation as animation
import matplotlib

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

#VOLTAMMETRY PARAMETERS
Cb0 = 1000 #bulk concentration of reduced species
sr = 0.008      #scan rate [V/s]
upl = 0.35     #upper potential limit [V]
sp = -0.35       #starting potential [V]
Ca0 = Cb0*np.exp(sp*frt) #concentration of A in equilibrium at starting potential


#ELECTRODE KINETIC PARAMETERS AND DIFFISION COEFFICIENT
C = 80          #double layer capacitance [F/m^2]
D = 4.0e-9      #Diffusion coefficient, m2/s

#TIME INTEGRATION PARAMETERS
t_start = 0                 #always start at time zero
#t_final = 4*(upl-sp)/sr     #completes just one cycle
t_final = 2*(upl-sp)/sr     #completes just one cycle
t_peak  = (upl-sp)/sr       #time to reach peak potential
t_min   = 2*(upl-sp)/sr     #time to reach minimum potential
dt = t_final/500.0          #step size for time
t = np.linspace(t_start, t_final, num=500)

pause = False

#ESTABLISH POTENTIAL PROFILE FOR SINGLE CV+ final anodic sweep
V = np.ones(len(t))
cd = np.zeros(len(t)) #set up array to record i(t), current density
s1 = int((len(t))/2)
s2 = int((len(t)))
for i in range (0, s1):
    V[i] = sp+sr*t[i]
for j in range(s1, s2):
    V[j] = upl-(t[j]-t_peak)*sr
for k in range(s2, len(t)):
    V[k] = -upl+(t[k]-t_min)*sr

"""Function to create charging current"""
idl = np.zeros(len(t))
for i in range (0, s1):
    idl[i] = sr*C
for j in range(s1, s2):
    idl[j] = -sr*C
for k in range(s2, len(t)):
    idl[k] = sr*C
    
#NUMERICAL PARAMETERS
N = 81          #number of grid points
Lc = 3*np.sqrt(2*D*t_final) #characteristic length, penetration depth for diffusion
nsteps = 499    #number of time steps
dx = Lc/(N-1)   #grid spacing
r = D*dt/dx**2
conc_A = np.zeros( (N, len(t)) )
conc_B = np.zeros( (N, len(t)) )

#INITIALIZE MATRICES A, B and b
A = np.zeros((N-2,N-2)) #python starts at zero, so we only go to N-2
B = np.zeros((N-2,N-2))
b = np.zeros((N-2))
#define matrices A, B and b array
for i in range(N-2):
    if i==0:
        A[i,:] = [2+2*r if j==0 else (-r) if j==1 else 0 for j in range(N-2)]
        B[i,:] = [2-2*r if j==0 else r if j==1 else 0 for j in range(N-2)]
        b[i] = 2*r*V[0] #boundary condition at i=1
    elif i==N-2:
        A[i,:] = [-r if j==N-4 else 2*2*r if j==N-3 else 0 for j in range(N-2)]
        B[i,:] = [r if j==N-4 else 2-2*r if j==N-3 else 0 for j in range(N-2)]
        b[i] = 0.0 #boundary condition at i=N
    else:
        A[i,:] = [-r if j==i-1 or j==i+1 else 2+2*r if j==i else 0 for j in range(N-2)]
        B[i,:] = [r if j==i-1 or j==i+1 else 2-2*r if j==i else 0 for j in range(N-2)]
#initialize grid
x = np.linspace(0,1,N)
#intial conditions
u = np.asarray([0.0 for xx in x])
w = np.asarray([0.0 for xx in x])
#u =  np.asarray([1.0 if xx<=0.5 else 2*(1-xx) for xx in x])

#evaluate right hand side at t=0
bb = B.dot(u[1:-1])+ b #calculate right hand side based on initial conditions
#syntax is start with 2nd element in matrix and stop at the next to the last one, the idea
#is that the first and last are determined by boundary conditions

matplotlib.rc('font', size=12)
matplotlib.rc('font', family='Arial')

#TIME MARCHING IMPLICIT SOLUTION
c = 0  #counter to keep track of step number
for j in range(nsteps): #j=0 corresponds to the first time step
#   print(j)
#find solution inside domain
    u[1:-1] = np.linalg.solve(A,bb)
#   u[0] =(1000.0 + 1000.0)*np.exp(frt*V[j])/(1.0 + np.exp(frt*V[j]))-1000.0
    u[0] =(Ca0 + Cb0)*np.exp(frt*V[j])/(1.0 + np.exp(frt*V[j]))
    for kk in range (N-1):
        w[kk]=-u[kk]
    
#update right hand side
    bb = B.dot(u[1:-1]) + b
    ca = (0.0012 + 1000.0)*np.exp(frt*V[j+1])/(1.0 + np.exp(frt*V[j+1]))
#   ca = (1000.0 + 1000.0)*np.exp(frt*V[j+1])/(1.0 + np.exp(frt*V[j+1]))-1000.0
    b[0] = 2*r*ca #boundary condition at i=1
#use fwd difference to calculate the current density at surface
    cd[j+1] = -F*D*(-u[2]+4*u[1]-3*u[0])/dx/2 + idl[j+1]
    conc_A[:,j]= Ca0 + u[:]
    conc_B[:,j]= Cb0 - u[:]
c +=1

#######################   ANIMATION
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8,8))
line1, = ax1.plot(t, V, 'r')
line2, = ax2.plot(V, cd, 'b') #CV
line3, = ax3.plot(t, cd, 'g')
line4, = ax4.plot(x, conc_A[:, 0], 'k', label=r'$c_A$')#concentration profile
line4b, = ax4.plot(x, conc_B[:, 0], 'b--', label=r'$c_B$')


def animate(i):
    if not pause:
        line1.set_xdata(t[1:i]) 
        line1.set_ydata(V[1:i])
        line2.set_ydata(cd[1:i])
        line2.set_xdata(V[1:i]) # update the data.
        line3.set_xdata(t[1:i]) 
        line3.set_ydata(cd[1:i])
        line4.set_ydata(conc_A[:, i])
        line4b.set_ydata(conc_B[:,i])
        return

def onClick(event):
    global pause
    pause ^= True
    
fig.canvas.mpl_connect('button_press_event', onClick)
ani = animation.FuncAnimation(fig, animate, frames=499, repeat=False, interval=3)
#ani = animation.FuncAnimation(fig, animate, frames=499, repeat=False, interval=320)

fig.suptitle('Classic Voltammogram, infinite kinetics', fontsize=18, y=0.93)
ax3.sharex(ax1)
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
ax4.grid(True)
ax2.set_xlim([-0.35, 0.35])
ax2.set_ylim([-2000,2000])
ax2.set_xticks(np.arange(-0.3, 0.35, step=0.1))
ax3.set_xlabel(r'Time, s')
ax4.set_xlabel(r'Distance from electrode')
ax4.set_ylim(0, 1500)
ax4.set_xlim(0,1)
plt.minorticks_on()
plt.locator_params(axis='x', nbins=10)
ax1.text(60, -0.18, 'Potential, V')
ax3.text(5, -400, 'Current density, ' r'$\mathrm {A m^{-2}}$')
ax4.text(0.5, 1200, 'A + ' r'$e^{-}$' ' =  B')
ax4.legend(bbox_to_anchor=(0.65, 0.6))

plt.show()
################################






























        



