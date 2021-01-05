"""
In this example, a slider is used to choose the concentration of sulfuric acid
and either a SHE or Hg-Hg2SO4 reference electrode can be selected.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.ticker import (MultipleLocator)

m0 = 1.0 #intial molality of acid

#Set up plot
fig, ax = plt.subplots(figsize=(8, 8))
plt.subplots_adjust(left=0.15, bottom=0.25)
ax.text(0.05, 2.6, 'Equilibrium potential for a lead-acid cell at 25$^\circ$C', fontsize=18)
ax.text(-0.1, 0.5, 'Potential, V', fontsize=18, rotation=90)
ax.text(0.28, 2.01, r'$U^\theta$', fontsize=12, color = 'k')
ax.text(-0.05, -0.365, r'$U^\theta$', fontsize=12, color = 'b')
ax.text(1.02, 1.65, r'$U^\theta$', fontsize=12, color = 'r')
ax.text(0.45, 0.7, r'full cell', fontsize=12)
ax.text(0.04, 0.7, r'negative electrode', fontsize=12, color='b')
ax.text(0.73, 0.7, r'positive electrode', fontsize=12, color = 'r')
ax.text(0.35, -1.8, r'molality of sulfuric acid, mol/kg', fontsize=12)
plt.xlim(0.0,1.0)
plt.ylim(-1.0,2.5)
ax.set_xticklabels([])
plt.grid(True)

#break up plot into three areas, negative, full, positive
ax.xaxis.set_major_locator(MultipleLocator(0.333))
t2 = np.arange(0.0, 0.33, 0.001)
t3 = np.arange(0.34, 0.66, 0.001)
t4 = np.arange(0.67, 1.0, 0.001)

#plot standard potentials
s2e = -0.356 + 0.0*t2 #standard potential of negative electrode
l2e, = plt.plot(t2, s2e, lw=3, color='b')#plot std potential of negative electrode
s3e = 2.041 + 0*t3 #standard potential of full cell
l3e, = plt.plot(t3, s3e, lw=3, color='k')
s4e = 1.685 + 0.0*t4#standard potential of positive electrode
l4e, = plt.plot(t4, s4e, lw=3, color='r')


#DEFINE FUNCTIONS for potentials and activity coefficients

def actcoef(m): #model for activity coefficient of acid
    q=np.sqrt(m)
    gam = 0.3368*q**6 -2.8004*q**5 + 9.138*q**4 - 14.783*q**3 + 12.38*q**2 - 5.1013*q+0.9633
    return gam

def actwat(m): #model for activity coefficient of water
    aw =1.0 - 0.0107*m -0.0023*m**2
    return aw

def pos(t4,m):
    gamm = actcoef(m)
    aw = actwat(m)
    s4 = 1.685 + 0.0385*np.log(m*gamm) - 0.0251*np.log(aw) + 0.0178 + 0.0*t4# potential of positive electrode
    return s4

def neg(t2,m):
    gamm = actcoef(m)
    s2 = -0.356 -0.0385*np.log(m*gamm) - 0.0178 + 0.0*t2# based on initial molality
    return s2

def ful(t3,m):
    gamm = actcoef(m)
    aw = actwat(m)
    s3 = 2.041 + 0.0771*np.log(m*gamm) - 0.0251*np.log(aw) + 0.0356 + 0.0*t3 #spotential of full cell
    return s3

#plot values for baseline conditions, 1 molal

s2 = neg(t2,m0)
l2, = plt.plot(t2, s2, lw=2, ls='--', color='b')

#s3 = 2.041 + 0.0771*np.log(m0*gamm) - 0.0251*np.log(aw) + 0.0356 + 0.0*t3 #spotential of full cell
s3 = ful(t3,m0)
l3, = plt.plot(t3, s3, lw=2, ls='--', color='k')

s4 = pos(t4,m0)
l4, = plt.plot(t4, s4, lw=2, ls='--', color='r')

#add cell potential text box
time_text = ax.text(0.45, 1.3, '',bbox=dict(facecolor='aliceblue'))
time_text.set_text(r'$U$ = %.3f'%  s3[2])

#radio button to select type of reference electrode
rax = plt.axes([0.45, 0.30, 0.16, 0.16], facecolor='aliceblue')
radio = RadioButtons(rax, ('SHE', 'Hg-Hg2SO4'))

axm = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor='aliceblue')
sm = Slider(axm, r'$m$', 0.001, 6.0, valinit=m0, color='dodgerblue', valfmt='%1.3f')

def update(val):#updates graph when slider moved
    tg = radio.value_selected
    m = sm.val
    gamm = actcoef(m)
    aw = actwat(m)
    #positive and negative electrode depend on selection of reference electrode    
    h2n = -0.356 -0.0385*np.log(m*gamm) - 0.0178 + 0.0*t2
    Hgn = -0.9717+ 0.0*t2+ 0.0*t2
    ren = {'SHE' : h2n, 'Hg-Hg2SO4' : Hgn}
    ydata = ren[tg]
    l2.set_ydata(ydata)
    h2p = 1.685 + 0.0385*np.log(m*gamm) - 0.0251*np.log(aw) + 0.0178 + 0.0*t4
    Hgp = 1.069 + 0.0771*np.log(m*gamm) - 0.0251*np.log(aw) + 0.0178 + 0.0*t4
    rep = {'SHE' : h2p, 'Hg-Hg2SO4' : Hgp}
    ydata = rep[tg]
    l4.set_ydata(ydata)
    #full cell is independent of reference electrode selection
    l3.set_ydata(2.041 + 0.0771*np.log(m*gamm) - 0.0251*np.log(aw) + 0.0356 + 0.0*t3)
    
    Ufs = 2.041 + 0.0771*np.log(m*gamm) - 0.0251*np.log(aw) + 0.0356
    time_text.set_text(r'$U$ = %.3f'% Ufs)
    fig.canvas.draw_idle()
    
def reff(label):#updates graph when reference electrode selected
    m = sm.val
    gamm = actcoef(m)
    aw = actwat(m)
    h2n = -0.356 -0.0385*np.log(m*gamm) - 0.0178 + 0.0*t2
    Hgn = -0.9717+ 0.0*t2+ 0.0*t2
    ren = {'SHE' : h2n, 'Hg-Hg2SO4' : Hgn}
    ydata = ren[label]
    l2.set_ydata(ydata)
    h2p = 1.685 + 0.0385*np.log(m*gamm) - 0.0251*np.log(aw) + 0.0178 + 0.0*t4
    Hgp = 1.069 + 0.0771*np.log(m*gamm) - 0.0251*np.log(aw) + 0.0178 + 0.0*t4
    rep = {'SHE' : h2p, 'Hg-Hg2SO4' : Hgp}
    ydata = rep[label]
    l4.set_ydata(ydata)
    fig.canvas.draw_idle()
    
sm.on_changed(update)

radio.on_clicked(reff)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color='silver', hovercolor='0.975')

def reset(event):#resets slider back to initial molality
    sm.reset()
button.on_clicked(reset)
   
plt.show()
