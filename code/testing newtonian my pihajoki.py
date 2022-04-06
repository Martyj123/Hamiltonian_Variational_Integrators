import math as ma
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
from scipy.constants import G
import time

#Define constants/initial values/arrays to store stuff in
m=100             #mass in kg
mm=100
hnx=10           #Initial position
hny=10
hnz=10
pnx=10           #Initial momentum 
pny=10
pnz=5
hnxt=hnx
hnyt=hny
hnzt=hnz
pnxt=pnx
pnyt=pny
pnzt=pnz
dt=0.001          #Time step
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
def dhdq(hi, hx, hy, hz): 
    r = np.sqrt(hx**2 + hy**2 + hz**2)
    answer = mm * hi  / r**3
    return answer

def dhdp(pi):
    answer = pi / m
    return answer
start = time.time()
hnxhalf = hnx + dhdp(pnx) * dt
hnyhalf = hny + dhdp(pny) * dt
hnzhalf = hnz + dhdp(pnz) * dt
pnxhalf = pnx - dhdq(hnx, hnx, hny, hnz) * dt
pnyhalf = pny - dhdq(hny, hnx, hny, hnz) * dt
pnzhalf = pnz - dhdq(hnz, hnx, hny, hnz) * dt
hnx1 = hnx + dhdp(pnxhalf) * dt * 2
hny1 = hny + dhdp(pnyhalf) * dt * 2
hnz1 = hnz + dhdp(pnzhalf) * dt * 2
pnx1 = pnx - dhdq(hnxhalf, hnxhalf, hnyhalf, hnzhalf) * dt * 2
pny1 = pny - dhdq(hnyhalf, hnxhalf, hnyhalf, hnzhalf) * dt * 2
pnz1 = pnz - dhdq(hnzhalf, hnxhalf, hnyhalf, hnzhalf) * dt * 2
hnx1 = hnxhalf + dhdp(pnx1) * dt
hny1 = hnyhalf + dhdp(pny1) * dt
hnz1 = hnzhalf + dhdp(pnz1) * dt
pnx1 = pnxhalf - dhdq(hnx1, hnx1, hny1, hnz1) * dt
pny1 = pnyhalf - dhdq(hny1, hnx1, hny1, hnz1) * dt
pnz1 = pnzhalf - dhdq(hnz1, hnx1, hny1, hnz1) * dt
hnx = hnx1
hny = hny1
hnz = hnz1
pnx = pnx1
pny = pny1
pnz = pnz1
t+=dt
tvals.append(t)
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 150:
    hnxhalf = hnx + dhdp(pnx) * dt
    hnyhalf = hny + dhdp(pny) * dt
    hnzhalf = hnz + dhdp(pnz) * dt
    pnxhalf = pnx - dhdq(hnx, hnx, hny, hnz) * dt
    pnyhalf = pny - dhdq(hny, hnx, hny, hnz) * dt
    pnzhalf = pnz - dhdq(hnz, hnx, hny, hnz) * dt
    hnx1 = hnx + dhdp(pnxhalf) * dt * 2
    hny1 = hny + dhdp(pnyhalf) * dt * 2
    hnz1 = hnz + dhdp(pnzhalf) * dt * 2
    pnx1 = pnx - dhdq(hnxhalf, hnxhalf, hnyhalf, hnzhalf) * dt * 2
    pny1 = pny - dhdq(hnyhalf, hnxhalf, hnyhalf, hnzhalf) * dt * 2
    pnz1 = pnz - dhdq(hnzhalf, hnxhalf, hnyhalf, hnzhalf) * dt * 2
    hnx1 = hnxhalf + dhdp(pnx1) * dt
    hny1 = hnyhalf + dhdp(pny1) * dt
    hnz1 = hnzhalf + dhdp(pnz1) * dt
    pnx1 = pnxhalf - dhdq(hnx1, hnx1, hny1, hnz1) * dt
    pny1 = pnyhalf - dhdq(hny1, hnx1, hny1, hnz1) * dt
    pnz1 = pnzhalf - dhdq(hnz1, hnx1, hny1, hnz1) * dt
    hnx = hnx1
    hny = hny1
    hnz = hnz1
    pnx = pnx1
    pny = pny1
    pnz = pnz1
    t+=dt
    tvals.append(t)
    xpos.append(hnx1)
    ypos.append(hny1)
    zpos.append(hnz1)
    xmom.append(pnx)
    ymom.append(pny)
    zmom.append(pnz)

end = time.time()
print(end - start)
ax = plt.axes(projection='3d')
# Data for three-dimensional scattered points
zdata = zpos
xdata = xpos
ydata = ypos
#ax.scatter3D(xdata, ydata, zdata, c=zdata);
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.plot3D(xdata, ydata, zdata, 'gray')
    
#plt.plot(tvals[:], evals[:])