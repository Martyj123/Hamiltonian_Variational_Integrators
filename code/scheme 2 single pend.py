import math as ma
import numpy as np
import matplotlib.pyplot as plt
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
qnminhalf=qn
pnminhalf=pn
qvals2=[]        #Store q values
pvals2=[]        #Store p values
tvals2=[]        #Store t values
posxvals=[]     #Store position in x
posyvals=[]
evals2=[]     #Store position in y
error2=[]
i=0
#Define derivatives to be used
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def dhdq(x):
    return -m * g * l * np.sin(x)

def dhdp(x):
    return x / (m * l ** 2)

initial=-m*g*l*(1-np.cos(qn))+(1/2) * (pn**2 / (m * l**2))
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 1000:
    #Page 1
    qnhalf = qnminhalf + dhdp(pn) * dt
    pnhalf = pnminhalf - dhdq(qn) * dt
    qn1 = qn + dhdp(pnhalf) * dt
    pn1 = pn - dhdq(qnhalf) * dt
    qn=qn1
    pn=pn1
    qnminhalf=qnhalf
    pnminhalf=pnhalf
    t+=dt;
    tvals2.append(t)
    posx=l * ma.sin(qn1)
    posy=l * ma.cos(qn1)
    posxvals.append(posx)
    posyvals.append(posy)
    pe=-m*g*l*(1-np.cos(qn1))
    ke= (1/2) * (pn1**2 /( m * l**2))
    E=pe+ke
    qvals2.append(ke)
    pvals2.append(pe)
    evals2.append(E)
    error2.append((E-initial)/initial)

average=sum(evals2)/(1000/dt)
print(average)
print(initial)
print(average-initial)
print(((average-initial)/initial)*100)
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
#                              interval=10, blit=True, init_func=init)
    
##Plotting the energy values
#plt.plot(tvals2, evals2, '-', color="black", markersize=1, label="Total Energy")
#plt.plot(tvals2, qvals2, '--', color="black", markersize=1, label="Kinetic Energy")
#plt.plot(tvals2, pvals2, '-.', color="black", markersize=1, label="Potential Energy")

#Plotting the relative error
plt.plot(tvals2, error2, '.', color="black", markersize=1, label="Relative error")

plt.xlabel('Time (s)')
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.ylabel('Energy')
plt.tick_params(axis="x", direction="in", length=5, width=1)
plt.tick_params(axis="y", direction="in", length=5, width=1)
plt.legend(loc='upper left', frameon=False)