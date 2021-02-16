from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt

############################################################

#STEP 1, Generating the grid    
M = 80
a = 0.0 #bottom left
b = 1.0 #bottom right
c = 0.0 #top left
d = 1.0 #top right

hx = (b - a)/M
N = int((d - c)/hx)#keeps sizes the same

#assumes that electrodes are centered
top_frac=0.30 #fraction of top covered by electrode
bot_frac=0.30 #fraction of bottom covered by electrode
t1 = int(M*(1-top_frac)/2)
t2 = int(M*(1+top_frac)/2)
b3 = int(N*(1-bot_frac)/2)
b4 = int(N*(1+bot_frac)/2)
ls1 = 0.5-top_frac/2#used to draw electrodes in graph
ls2 = 0.5-bot_frac/2

Wa = 1.0
pm = 2*hx/(top_frac*Wa)

x1 = np.linspace(a, b, M+1)
y1 = np.linspace(c, d, N+1)

X, Y = np.meshgrid(x1, y1)#Make 2-D array of grid 

#STEP 2, Loading the grid
#----- Right hand side of differential equation 

# a. diagonal and off diagonal elements
main_diag = 4*np.ones((M+1, 1)).ravel() #alpha is equal to 1, so this is 4
off_diag = -1*np.ones((M, 1)).ravel()

a = main_diag.shape[0]

diagonals = [main_diag, off_diag, off_diag]

B = sparse.diags(diagonals, [0,-1,1], shape=(a,a)).toarray()
B[0,1] = -2 #changes one value from -1 to -2
B[M,M-1] = -2 #not in original code, mistake?

D = sparse.diags([-1*np.ones((M+1, 1)).ravel()], [0], shape=(a,a)).toarray()

C = sparse.diags([-2*np.ones((M+1, 1)).ravel()], [0], shape=(a,a)).toarray()

e1 = sparse.eye(M+1).toarray() # 2d array where diags is one

A1 = sparse.kron(e1,B).toarray() #here is where B is used

e2 = sparse.diags([1*np.ones((M, 1)).ravel(),1*np.ones((M, 1)).ravel()], [-1,1], shape=(M+1,M+1)).toarray()
e2[0,:] = 0
e2[M,:] = 0

A2 = sparse.kron(e2,D).toarray() #here is where D is used

e3 = sparse.diags([1*np.ones((M, 1)).ravel(),1*np.ones((M, 1)).ravel()], [-1,1], shape=(M+1,M+1)).toarray()
e3[1:M,0:M+1] = 0

A3 = sparse.kron(e3,C).toarray() #here is where C is used

A = A1 + A2 + A3

# The zero-flux boundary conditions are incorporated in matrix A.

g = np.zeros(((M+1)**2, 1))
#g[5:10] = g[5:10] + 5.0*hx


#convert part of bottom to specifying value, zero
bh = (M+1)*(M+1)-t1
bl = (M+1)*(M+1)-(t2)
for nb in range(bl, bh):
    A[nb,nb] = A[nb,nb] + pm
#   A[nb,:] = 0.0
#   A[nb,nb] = 1.0
#g[bl:bh] = 0.0
    
#convert part of bottom to specifying value, one
th = (M+1)-b3
tl = (M+1)-b4
for nb in range(tl, th): #with range, stops at th-1
    A[nb,:] = 0.0
    A[nb,nb] = 1.0
g[tl:th] = 0.10

U = np.linalg.solve(A,g)  # Solve x=A\b
U = U.reshape((M+1,M+1)).T #reshape back into 2d matrix


#PLOT OF ISOPOENTIAL AND STREAMLINES

#Calculate gradient of field
gradx, grady = np.gradient (U, hx, hx)

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(16,8))

ax1.set_xlim(0.0,1.0)
ax1.set_ylim(0.0,1.0)

ax1.streamplot(y1, x1, -grady, -gradx, color ="orange", linewidth=0.5, density=[0.6,0.6])
CS = ax1.contour(X, Y, U, [0.01, 0.04, 0.045, 0.05, 0.06, 0.07, 0.08, 0.085, 0.09], colors='b', linestyles='--') 
ax1.clabel(CS, fontsize=12, inline=0, fmt='%1.3f')
ax1.add_patch(plt.Rectangle((1.0,ls1),0.01, top_frac, facecolor='k',
                              clip_on=False,linewidth = 0))
ax1.add_patch(plt.Rectangle((-0.01,ls2),0.01, bot_frac,facecolor='k',
                              clip_on=False,linewidth = 0))

ax1.spines['left'].set_linewidth(1)
ax1.spines['right'].set_linewidth(1)
ax1.spines['top'].set_linewidth(1)
ax1.spines['bottom'].set_linewidth(1)
ax1.text(1.01, 0.50, 'Wa=', fontsize=14)
ax1.text(1.095, 0.50, Wa, fontsize=14)
ax1.text(-0.15, 0.50, 'Wa=0', fontsize=14)
ax2.text(0.22,0.28, 'Current density on the electrodes', fontsize=16)
ax1.tick_params(labelsize=16)

ax2.set_xlim(0.0,1.0)
#ax2.set_ylim(0.33, 0.67)
"""ax2.set_xlabel(r'Normalized variation in current density on electrode', fontsize=14)
ax1.set_xlabel(r'Dimensionless distance from left electrode', fontsize=14)
ax2.spines['left'].set_linewidth(1)
ax2.spines['right'].set_linewidth(1)
ax2.spines['top'].set_linewidth(1)
ax2.spines['bottom'].set_linewidth(1)
ax2.grid(True, which='both')"""

ax2.xaxis.set_visible(False)
ax2.yaxis.set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.spines["left"].set_visible(False)
ax2.spines["bottom"].set_visible(False)


xx = X[M,:]
yy = 1+grady[:,M]
yz = -grady[:,0]
ay = np.mean(yy[t1+1:t2+1])
az = np.mean(yz[t1+1:t2+1])
yyn = yy/ay/3 + 0.6 
yzn = yz/az/3 - 0.2
line1, = ax2.plot(yyn[t1+1:t2+1], xx[t1+1:t2+1], 'b')
line2, = ax2.plot(yzn[t1+1:t2+1], xx[t1+1:t2+1], 'r')



"""Wa_right_mean = np.mean(yyn[t1+1:t2+1])
Wa_left_mean  = np.mean(yzn[t1+1:t2+1])
Wa_right_variance = np.var(yyn[t1+1:t2+1])
Wa_left_variance  = np.var(yzn[t1+1:t2+1])
#print(Wa_left_mean, Wa_right_mean)"""
#print(Wa_left_variance, Wa_right_variance)
fig.suptitle('Secondary Current Distribution', fontsize=18, y=0.93)
anno1 = ax2.annotate('', xytext=(0,0.36), xy=(1.0,0.36), arrowprops=dict(color='k', lw=2, arrowstyle="-"))
anno2 = ax2.annotate('', xytext=(0,0.655), xy=(1.0,0.655), arrowprops=dict(color='k', lw=2, arrowstyle="-"))
anno3 = ax2.annotate('', xytext=(0,0.355), xy=(0.0,0.655), arrowprops=dict(color='k', lw=2, arrowstyle="-"))
anno4 = ax2.annotate('', xytext=(1.0,0.355), xy=(1.0,0.655), arrowprops=dict(color='k', lw=2, arrowstyle="-"))
plt.show()


