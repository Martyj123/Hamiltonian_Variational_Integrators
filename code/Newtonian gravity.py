import math as ma
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
from scipy.constants import G

#Define constants/initial values/arrays to store stuff in
m=1             #mass in kg
mm=1
hnx=1           #Initial position
hny=1
hnz=1
pnx=1           #Initial momentum 
pny=0
pnz=0
positionatn = [hnx, hny, hnz]
momentumatn = [pnx, pny, pnz]
positionathalf = positionatn
momentumathalf = momentumatn
positionatn1 = [0,0,0]
momentumatn1 = [0,0,0]
positionat3o2 = [0,0,0]
momentumat3o2 = [0,0,0]
dt=0.1          #Time step
t=0             #Initial time
tvals=[]        #Store t values
xpos=[]
ypos=[]
zpos=[]
xmom=[]
ymom=[]
zmom=[]
evals=[]
#Define derivatives to be used
def dhdq(x, i): 
    return mm * m  * i / (x[0]**2 + x[1]**2 + x[2]**2)**(3/2)

def dhdp(x):
    return x / m

#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 10:
    i=0
    while i<3:
        positionatn1[i] = positionatn[i] + dhdp(momentumathalf[i]) * dt
        momentumatn1[i] = momentumatn[i] - dhdq(positionathalf, positionathalf[i]) * dt
        positionat3o2[i] = positionathalf[i] + dhdp(momentumatn1[i]) * dt
        momentumat3o2[i] = momentumathalf[i] - dhdq(positionatn1, positionatn1[i]) * dt
        i+=1
    i=0
    while i<3:
        positionathalf[i] = positionat3o2[i]
        momentumathalf[i] = momentumat3o2[i]
        positionatn[i] = positionatn1[i]
        momentumatn[i] = momentumatn1[i]
        i+=1
    xpos.append(positionatn[0])
    ypos.append(positionatn[1])
    zpos.append(positionatn[2])
    xmom.append(momentumatn[0])
    ymom.append(momentumatn[1])
    zmom.append(momentumatn[2])
    t+=dt
    tvals.append(t)


ax = plt.axes(projection='3d')
# Data for three-dimensional scattered points
zdata = zpos
xdata = xpos
ydata = ypos
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens');
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
    
#plt.plot(tvals[:], evals[:])