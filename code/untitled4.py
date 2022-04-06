import math as ma
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

#Define constants/initial values/arrays to store stuff in
m=1             #mass in kg
l=10            #length in m
g=9.81         #Grav. acceleration in m/ss
qn=np.pi/4        #Initial angle of the pendulum
pn=0            #Initial momentum of the pendulum
dt=0.001         #Time step
t=0             #Initial time
qnm1=qn
pnm1=pn
qvals1=[]        #Store q values
pvals1=[]        #Store p values
tvals1=[]        #Store t values
posxvals=[]     #Store position in x
posyvals=[]
evals1=[]     #Store position in y
error1=[]
uncert=[]
ip=[]
i=20
#Define derivatives to be used
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def dhdq(x):
    return -m * g * l * ma.sin(x)

def dhdp(x):
    return x / (m * l ** 2)

initial=-m*g*l*(1-ma.cos(qn))+(1/2) * (pn**2 / (m * l**2))
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t<1000:
    qn=np.pi/4        #Initial angle of the pendulum
    pn=0            #Initial momentum of the pendulum
    dt=0.001         #Time step
    t=0             #Initial time
    qnm1=qn
    pnm1=pn
    qvals1=[]        #Store q values
    pvals1=[]        #Store p values
    tvals1=[]        #Store t values
    posxvals=[]     #Store position in x
    posyvals=[]
    evals1=[]     #Store position in y
    error1=[]
    while t < i:
        #Page 1
        qn1 = qnm1 + dhdp(pn) * dt * 2
        pn1 = pnm1 - dhdq(qn) * dt * 2
        qn2 = qn + dhdp(pn1) * dt * 2
        pn2 = pn - dhdq(qn1) * dt * 2
        qnm1=qn1
        pnm1=pn1
        qn=qn2
        pn=pn2
        posx=l * ma.sin(qn)
        posy=l * ma.cos(qn)
        posxvals.append(posx)
        posyvals.append(posy)
        pe=-m*g*l*(1-ma.cos(qn))
        ke= (1/2) * (pn**2 / (m * l**2))
        E=pe+ke
        qvals1.append(qn)
        pvals1.append(pn)
        evals1.append(E)
        error1.append((E-initial)/initial)
        t+=dt
        tvals1.append(t)
    
    average=sum(evals1)/(t/dt)
    uncertainty=((average-initial)/initial)*100
    uncert.append(uncertainty)
    print(abs(uncertainty))
    ip.append(i)
    i+=10
    
#Plotting the energy values
#plt.plot(tvals1, evals1, '.', color="black", markersize=1, label="Total Energy")
#plt.plot(tvals1, qvals1, '--', color="black", markersize=1, label="Kinetic Energy")
#plt.plot(tvals1, pvals1, '-.', color="black", markersize=1, label="Potential Energy")

#Plotting the relative error
plt.plot(ip, uncert, '-', color="black", linewidth=1, label="Relative error")
plt.xlabel('Time (s)')
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.ylabel('Relative error')
plt.tick_params(axis="x", direction="in", length=5, width=1)
plt.tick_params(axis="y", direction="in", length=5, width=1)
plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 3))
#plt.xscale("log")
#plt.yscale("log")

#fig = plt.figure()
#ax = fig.add_subplot(111, autoscale_on=False, xlim=(-15, 15), ylim=(-15, 15))
#ax.set_aspect('equal')
#ax.grid()
#
#line, = ax.plot([], [], 'o-', lw=2)
#time_template = 'time = %.1fs'
#time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
#
#
#def animate(i):
#    thisx = [0, posxvals[i]]
#    thisy = [0, posyvals[i]]
#    line.set_data(thisx, thisy)
#    time_text.set_text(time_template % (i*dt))
#    return line, time_text
#
#ani = animation.FuncAnimation(fig, animate,
#                              interval=1, blit=True, init_func=init)