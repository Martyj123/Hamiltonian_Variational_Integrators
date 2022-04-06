import math as ma
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

#Define constants/initial values/arrays to store stuff in
m=1             #mass in kg
l=10            #length in m
g=9.81         #Grav. acceleration in m/ss
qn=-3        #Initial angle of the pendulum
pn=0            #Initial momentum of the pendulum
dt=0.001         #Time step
t=0             #Initial time
qnm1=qn
pnm1=pn
qvals1=[]        #Store q values
pvals1=[]        #Store p values
tvals2v=[]        #Store t values
posxvals=[]     #Store position in x
posyvals=[]
evals2v=[]     #Store position in y
error2v=[]
i=0
#Define derivatives to be used
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def hamilt(x,y):
    return (1/2)*((x**2)*(y**2) +x**2+y**2+1)

def dhdq(x,y):
    return x*y**2 + x

def dhdp(x,y):
    return y*x**2 + y

initial=hamilt(qn,pn)
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 10:
    #Page 1
    q1half = qn + dhdp(qn, pn) * dt /2
    p1half = pn - dhdq(qn, pn) * dt /2
    q1n1 = qn + dhdp(q1half, p1half) * dt 
    p1n1 = pn - dhdq(q1half, p1half) * dt 
    q1n1 = q1half + dhdp(q1n1, p1n1) * dt /2
    p1n1 = p1half - dhdq(q1n1, p1n1) * dt /2
    #Update values for next time step
    qn = q1n1
    pn = p1n1
    t+=dt
    tvals2v.append(t)
    posx=l * ma.sin(qn)
    posy=l * ma.cos(qn)
    posxvals.append(posx)
    posyvals.append(posy)
    E=hamilt(qn,pn)
    qvals1.append(qn)
    pvals1.append(pn)
    evals2v.append(E)
    error2v.append((E-initial)/initial)


average=sum(evals2v)/(10/dt)
print(average)
print(initial)
print(average-initial)
print(((average-initial)/initial)*100)

#Plotting the energy values
#plt.plot(tvals1, evals1, '.', color="black", markersize=1, label="Total Energy")
#plt.plot(tvals1, qvals1, '--', color="black", markersize=1, label="Kinetic Energy")
#plt.plot(tvals1, pvals1, '-.', color="black", markersize=1, label="Potential Energy")

#plt.subplot(1,3,1)
#plt.title("Scheme 1")
#plt.ylabel('Relative error')
#plt.plot(tvals1n, error1n, '-', color="black", linewidth=1, label="Relative error")
#plt.xlabel('Time (s)')
#plt.tick_params(bottom=True, top=True, left=True, right=True)
#plt.tick_params(axis="x", direction="in", length=5, width=1)
#plt.tick_params(axis="y", direction="in", length=5, width=1)
#plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 3))
#
#plt.subplot(1,3,2)
#plt.title("Scheme 2")
#plt.plot(tvals2n, error2n, '-', color="black", linewidth=1, label="Relative error")
#plt.xlabel('Time (s)')
#plt.tick_params(bottom=True, top=True, left=True, right=True)
#plt.tick_params(axis="x", direction="in", length=5, width=1)
#plt.tick_params(axis="y", direction="in", length=5, width=1)
#plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 3))
#
#plt.subplot(1,3,3)
plt.title("PSVM")
plt.plot(tvals2v, error2v, '-', color="black", linewidth=1, label="Relative error")
plt.xlabel('Time (s)')
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.tick_params(axis="x", direction="in", length=5, width=1)
plt.tick_params(axis="y", direction="in", length=5, width=1)
plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 3))