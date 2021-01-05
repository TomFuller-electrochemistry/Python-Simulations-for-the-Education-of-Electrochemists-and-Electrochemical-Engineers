import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button
import matplotlib
#import scipy.integrate as integrate
"""import csv
fields = ['time', 'potential', 'current density'] 
filename = "university_records.csv" """
#SET UP PLOT
matplotlib.rc('font', size=14)
plt.rc('font',family='Times New Roman')
plt.rc('font', family='Arial')
fig1, ax = plt.subplots(figsize=(7, 7))
plt.xlabel('Overpotential, V', size=18)
plt.ylabel('Current density, A-$\mathrm{m^{-2}}$', size=18)
plt.text(-0.50, 650, 'A + ' r'$e^{-}$' ' =  B')
plt.text(-0.50, 500, r'CV for inifinite kinetics')
plt.xlim([-0.6, 0.6])
plt.ylim([-1000,1000])
plt.grid(True)
plt.minorticks_on()
plt.locator_params(axis='x', nbins=10)
plt.subplots_adjust(left=0.15)

#TEMPERATURE AND CONSTANTS
T_ex = 298.15   #Temperature, K
F = 96485.0     #Faraday's constant
R = 8.414       #Universal gas constant
frt = F/R/T_ex

#VOLTAMMETRY PARAMETERS
sr0 = 0.005    #scan rate [V/s]
sr = sr0
upl = 0.50  #upper potential limit [V]

#ELECTRODE KINETIC PARAMETERS AND DIFFISION COEFFICIENT
C0 = 2        #double layer capacitance [F/m^2]
C = C0
D0 = 1.0e-9      #Diffusion coefficient, m2/s
D = D0
cr = 500

#Slider for three parameters: scan rate, capacitance, and diffusivity
axm = plt.axes([0.35, 0.885, 0.40, 0.02], facecolor='aliceblue')
srm = Slider(axm, r'scan rate [mV/s]', 1, 20, valinit=1000*sr0, color='dodgerblue', valfmt='%1.1f')
axm2 = plt.axes([0.35, 0.915, 0.40, 0.02], facecolor='aliceblue')
cap = Slider(axm2, r'capacitance $[\mathrm{F/ cm^2]}$', 0, 4, valinit=C0/1e4, color='dodgerblue', valfmt='%1.1f')
axm3 = plt.axes([0.35, 0.945, 0.40, 0.02], facecolor='aliceblue')
diff = Slider(axm3, r'diffusivity x $10^9 [\mathrm{m^2/s}]$', 0.1, 10, valinit=D0*1e9, color='dodgerblue', valfmt='%1.1f')
#srm.label.set_size(12)

#DEFINE FUNCTIONS
def volt(sr, t, t_peak, t_min):#sets up voltage vs. time based on scane rate
    for i in range (0, s1):
        V[i] = 0.0 + sr*t[i]
    for j in range(s1, s2):
        V[j] = upl-(t[j]-t_peak)*sr
    for k in range(s2, len(t)):
        V[k] = -upl+(t[k]-t_min)*sr
    return V
          
def intpara(scr): #TIME INTEGRATION PARAMETERS
    D = diff.val/1e9
    t_final = 4*(upl)/scr     #completes one cycle
    t_peak  = (upl)/scr       #time to reach peak potential
    t_min   = 3*(upl)/scr     #time to reach minimum potential
    dt = t_final/500.0         #step size for time
    Lc = 3*np.sqrt(2*D*t_final) #characteristic length, penetration depth for diffusion
    dx = Lc/(N-1)   #grid spacing
    r = D*dt/dx**2
    return t_final, t_peak, t_min, dt, Lc, dx, r
         
def chcur(sr, C):
    idl = np.zeros(len(t))
    for i in range (0, s1):
        idl[i] = sr*C
    for j in range(s1, s2):
        idl[j] = -sr*C
    for k in range(s2, len(t)):
        idl[k] = sr*C
    return idl

def update(val):#updates graph when slider moved
    sr = srm.val/1000
    C = cap.val*1e4
    D = diff.val/1e9
    t_final, t_peak, t_min, dt, Lc, dx, r = intpara(sr)
    t = np.linspace(0, t_final, num=500)
    V = volt(sr, t, t_peak, t_min)
    idl = chcur(sr, C)
    scan(sr, bb, idl, V, ca, A, B, b, dx, D)
    x = V
    y = cd
    line.set_xdata(x[1:])
    line.set_ydata(y[1:])
    fig1.canvas.draw_idle()
"""    with open(filename, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        rows = [t[1:], V[1:], y[1:]]
        csvwriter.writerows(rows)"""

#    coulpass = integrate.simps(yint, tint)
#   line, = ax.plot(x, y, 'g-o')



def scan(sr, bb, idl, V, ca, A, B, b, dx, D):
    c = 0  #counter to keep track of step number
    for j in range(nsteps): #j=0 corresponds to the first time step
        #   print(j)
        #find solution inside domain
        u[1:-1] = np.linalg.solve(A,bb)
        u[0] =(cr + cr)*np.exp(frt*V[j])/(1.0 + np.exp(frt*V[j]))-cr#concentration
        for kk in range (N-1):
            w[kk]=-u[kk]
            
        #update right hand side
        bb = B.dot(u[1:-1]) + b
        ca = (cr + cr)*np.exp(frt*V[j+1])/(1.0 + np.exp(frt*V[j+1]))-cr#concentration
        b[0] = 2*r*ca #boundary condition at i=1
        #use fwd difference to calculate the current density at surface
        cd[j+1] = -F*D*(-u[2]+4*u[1]-3*u[0])/dx/2 + idl[j+1]
        
    c +=1

#   line, = ax.plot(x, y, 'g-o')

def animate(i):
    line.set_ydata(cd[1:i])
    line.set_xdata(V[1:i]) # update the data.
#    if i > 499:
#         import sys
#         sys.exit()
#    line.set = ax.plot(x, y, 'g-o')
    return line,

#INTEGRATION and NUMERICAL PARAMETERS based on original scan rate
t_final = 4*(upl)/sr     #completes one cycle
t_peak  = (upl)/sr       #time to reach peak potential
t_min   = 3*(upl)/sr     #time to reach minimum potential
dt = t_final/500.0          #step size for time

N = 121          #number of grid points
Lc = 4*np.sqrt(2*D*t_final) #characteristic length, penetration depth for diffusion
nsteps = 499    #number of time steps
dx = Lc/(N-1)   #grid spacing
#bc = 0.0        #boundary condition at x=0
r = D*dt/dx**2

t_final, t_peak, t_min, dt, Lc, dx, r = intpara(sr) 

t = np.linspace(0, t_final, num=500)

#ESTABLISH POTENTIAL PROFILE FOR SINGLE CV+ final anodic sweep
V = np.ones(len(t))
cd = np.zeros(len(t)) #set up array to record i(t), current density
s1 = int((len(t))/4)
s2 = int((3*len(t))/4)

V = volt(sr, t, t_peak, t_min)#set potential as a function of time

"""Calculate charging current"""
idl = np.zeros(len(t))
idl = chcur(sr, C)
    

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

#TIME MARCHING SOLUTION
c = 0  #counter to keep track of step number
for j in range(nsteps): #j=0 corresponds to the first time step
#   print(j)
#find solution inside domain
    u[1:-1] = np.linalg.solve(A,bb)
    u[0] =(cr + cr)*np.exp(frt*V[j])/(1.0 + np.exp(frt*V[j]))-cr#concentration
    for kk in range (N-1):
        w[kk]=-u[kk]
    
#update right hand side
    bb = B.dot(u[1:-1]) + b
    ca = (cr + cr)*np.exp(frt*V[j+1])/(1.0 + np.exp(frt*V[j+1]))-cr
    b[0] = 2*r*ca #boundary condition at i=1
#use fwd difference to calculate the current density at surface
#   print(D, r, dx)
    cd[j+1] = -F*D*(-u[2]+4*u[1]-3*u[0])/dx/2 + idl[j+1]
c +=1

x = V
y = cd
#line, = ax.plot(x, y, 'r-o')
fstscan, = ax.plot(x, y, 'r')

#ani = animation.FuncAnimation(fig1, animate, interval=10)
#plt.show()

#SIMPLY REPEATING THE INTEGRATION FOR SECOND CYCLE
c = 0  #counter to keep track of step number
for j in range(nsteps): #j=0 corresponds to the first time step
#   print(j)
#find solution inside domain
    u[1:-1] = np.linalg.solve(A,bb)
    u[0] =(cr + cr)*np.exp(frt*V[j])/(1.0 + np.exp(frt*V[j]))-cr#concentration
    for kk in range (N-1):
        w[kk]=-u[kk]
    
#update right hand side
    bb = B.dot(u[1:-1]) + b
    ca = (cr+ cr)*np.exp(frt*V[j+1])/(1.0 + np.exp(frt*V[j+1]))-cr#cr
    b[0] = 2*r*ca #boundary condition at i=1
#use fwd difference to calculate the current density at surface
    cd[j+1] = -F*D*(-u[2]+4*u[1]-3*u[0])/dx/2 + idl[j+1]

c +=1

x = V
y = cd

line, = ax.plot(x[1:], y[1:], lw=2, color='k')

ani = animation.FuncAnimation(fig1, animate, interval=20)

#additional cycles initiated on movement of sliders
srm.on_changed(update)
cap.on_changed(update)
diff.on_changed(update)

#anmax = plt.axes([0.65, 0.30, 0.15, 0.05])
#button1 = Button(anmax, 'Animate', color='silver', hovercolor='0.975')

fstx = plt.axes([0.65, 0.20, 0.22, 0.05])
button2 = Button(fstx, 'Remove 1st scan', color='red', hovercolor='0.975')

"""def dynscan(event):
#    x = V
#    y = cd
    line, = ax.plot(V, cd, 'w-o')
    ani = animation.FuncAnimation(fig1, animate, interval=20, blit=True)
    return line,

#   button1.set_active(False)"""


def removefirstscan(event):
    fstscan.remove()
#   button2.set_active(False)
button2.on_clicked(removefirstscan)    
#button1.on_clicked(dynscan)


################################































        



