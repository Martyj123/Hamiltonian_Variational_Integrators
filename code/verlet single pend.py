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
qvalsv=[]        #Store q values
pvalsv=[]        #Store p values
tvalsv=[]        #Store t values
posxvals=[]     #Store position in x
posyvals=[]
evalsv=[]     #Store position in y
errorv=[]
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
while t < 100:
    #Page 1
    qnhalf = qn + dhdp(pn) * dt /2 
    pnhalf = pn - dhdq(qn) * dt/2
    pn1 = pn - dhdq(qnhalf) * dt
    qn1 = qn + dhdp(pnhalf) * dt
    pn1n = pnhalf - dhdq(qn1) * dt/2
    qn1n = qnhalf + dhdp(pn1) * dt/2
    qn=qn1n
    pn=pn1n
    t+=dt;
    tvalsv.append(t)
    posx=l * ma.sin(qn1n)
    posy=l * ma.cos(qn1n)
    posxvals.append(posx)
    posyvals.append(posy)
    pe=-m*g*l*(1-np.cos(qn1n))
    ke= (1/2) * (pn1n**2 / (m * l**2))
    E=pe+ke
    qvalsv.append(ke)
    pvalsv.append(pe)
    evalsv.append(E)
    errorv.append((E-initial)/initial)

average=sum(evalsv)/(100/dt)
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
#                              interval=1, blit=True, init_func=init)


#plt.subplot(1, 3, 1)
#plt.title("Scheme 1")
###Plotting the energy values
#plt.plot(tvals1, evals1, '-', color="black", markersize=1, label="Total")
#plt.plot(tvals1, qvals1, '--', color="red", markersize=1, label="Kinetic")
#plt.plot(tvals1, pvals1, ':', color="blue", markersize=1, label="Potential")
#plt.xlabel('Time (s)')
#plt.tick_params(bottom=True, top=True, left=True, right=True)
#plt.ylabel('Energy (J)')
#plt.tick_params(axis="x", direction="in", length=5, width=1)
#plt.tick_params(axis="y", direction="in", length=5, width=1)
#
#plt.subplot(1, 3, 2)
#plt.title("Scheme 2")
###Plotting the energy values
#plt.plot(tvals2, evals2, '-', color="black", markersize=1, label="Total")
#plt.plot(tvals2, qvals2, '--', color="red", markersize=1, label="Kinetic")
#plt.plot(tvals2, pvals2, ':', color="blue", markersize=1, label="Potential")
#plt.xlabel('Time (s)')
#plt.tick_params(bottom=True, top=True, left=True, right=True)
#plt.tick_params(axis="x", direction="in", length=5, width=1)
#plt.tick_params(axis="y", direction="in", length=5, width=1)
#
#plt.subplot(1, 3, 3)
#plt.title("Stormer-Verlet scheme")
###Plotting the energy values
#plt.plot(tvalsv, evalsv, '-', color="black", markersize=1, label="Total")
#plt.plot(tvalsv, qvalsv, '--', color="red", markersize=1, label="Kinetic")
#plt.plot(tvalsv, pvalsv, ':', color="blue", markersize=1, label="Potential")
#plt.xlabel('Time (s)')
#plt.tick_params(bottom=True, top=True, left=True, right=True)
#plt.tick_params(axis="x", direction="in", length=5, width=1)
#plt.tick_params(axis="y", direction="in", length=5, width=1)
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
#plt.subplot(1,3,1)
#plt.title("Scheme 1")
#plt.ylabel('Relative error')
#plt.plot(tvals1, error1, '-', color="black", linewidth=1, label="Relative error")
#plt.xlabel('Time (s)')
#plt.tick_params(bottom=True, top=True, left=True, right=True)
#plt.tick_params(axis="x", direction="in", length=5, width=1)
#plt.tick_params(axis="y", direction="in", length=5, width=1)
#plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 3))
#
#plt.subplot(1,3,2)
#plt.title("Scheme 2")
#plt.plot(tvals2, error2, '-', color="black", linewidth=1, label="Relative error")
#plt.xlabel('Time (s)')
#plt.tick_params(bottom=True, top=True, left=True, right=True)
#plt.tick_params(axis="x", direction="in", length=5, width=1)
#plt.tick_params(axis="y", direction="in", length=5, width=1)
#plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 3))
#
#plt.subplot(1,3,3)
plt.title("PSVM")
plt.plot(tvalsv, errorv, '-', color="black", linewidth=1, label="Relative error")
plt.xlabel('Time (s)')
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.tick_params(axis="x", direction="in", length=5, width=1)
plt.tick_params(axis="y", direction="in", length=5, width=1)
plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 3))