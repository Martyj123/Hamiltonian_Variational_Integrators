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
pn=0.1            #Initial momentum of the pendulum
dt=0.001         #Time step
t=0             #Initial time
qnm1=qn
pnm1=pn
qvals1=[]        #Store q values
pvals1=[]        #Store p values
tvals2n=[]        #Store t values
posxvals=[]     #Store position in x
posyvals=[]
evals2n=[]     #Store position in y
error2n=[]
i=0
#Define derivatives to be used
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def hamilt(x,y):
    return (1/2)*((x**2)*(y**2) + (x**2) + (y**2) +1)

def dhdq(x,y):
    return x*y**2 + x

def dhdp(x,y):
    return y*x**2 + y

initial=hamilt(qn,pn)
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 100:
    #Page 1
    qn1 = qnm1 + dhdp(qn,pn) * dt
    pn1 = pnm1 - dhdq(qn,pn) * dt
    qn2 = qn + dhdp(qn1,pn1) * dt
    pn2 = pn - dhdq(qn1,pn1) * dt
    qnm1=qn1
    pnm1=pn1
    qn=qn2
    pn=pn2
    t+=dt
    tvals2n.append(t)
    E=hamilt(qn,pn)
    qvals1.append(qn)
    pvals1.append(pn)
    evals2n.append(E)
    error2n.append((E-initial)/initial)

maxerror=max(error2n)
average=sum(evals2n)/(100/dt)
print(average)
print(initial)
print(average-initial)
print(maxerror)
print(((average-initial)/initial)*100)

#Plotting the energy values
#plt.plot(tvals1, evals1, '.', color="black", markersize=1, label="Total Energy")
#plt.plot(tvals1, qvals1, '--', color="black", markersize=1, label="Kinetic Energy")
#plt.plot(tvals1, pvals1, '-.', color="black", markersize=1, label="Potential Energy")

#Plotting the relative error
plt.plot(tvals2n, error2n, '-', color="black", linewidth=1, label="Relative error")

plt.xlabel('Time (s)')
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.ylabel('Energy')
plt.tick_params(axis="x", direction="in", length=5, width=1)
plt.tick_params(axis="y", direction="in", length=5, width=1)
plt.legend(loc='upper left', frameon=False)
plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 3))