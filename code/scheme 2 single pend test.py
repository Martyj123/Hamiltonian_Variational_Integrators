import math as ma
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

#Define constants/initial values/arrays to store stuff in
m=1             #mass in kg
l=10            #length in m
g=-9.81         #Grav. acceleration in m/ss
qnminhalf=np.pi/4        #Initial angle of the pendulum
pnminhalf=0            #Initial momentum of the pendulum
dt=0.001         #Time step
t=0             #Initial time
qnminone = qnminhalf + qnminhalf * dt
qvals=[]        #Store q values
pvals=[]        #Store p values
tvals=[]        #Store t values
posxvals=[]     #Store position in x
posyvals=[]
evals=[]     #Store position in y
i=0
#Define derivatives to be used
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def dhdq(x):
    return m * g * l * np.sin(x)

def dhdp(x):
    return x / (m * l ** 2)
qn = qnminone + dhdp(pnminhalf) * dt
pn = -dhdq(qnminhalf) * dt
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 100:
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
    tvals.append(t)
    posx=l * ma.sin(qn1)
    posy=l * ma.cos(qn1)
    posxvals.append(posx)
    posyvals.append(posy)
    pe=m*g*l*(1-np.cos(qn1))
    ke= (1/2) * (pn1**2 /( m * l**2))
    E=pe+ke
    qvals.append(ke)
    pvals.append(pe)
    evals.append(E)

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
plt.xlabel('Time (s)')
plt.ylabel('Energy')
plt.plot(tvals[:], evals[:])
plt.plot(tvals[:], qvals[:])
plt.plot(tvals[:], pvals[:])