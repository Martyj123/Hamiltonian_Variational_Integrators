import math as ma
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
from scipy.constants import G
from scipy.constants import c
import time
#Define constants/initial values/arrays to store stuff in
m=1.659e-7           #mass in kg
mm=1
hnx=-0.46           #Initial position
hny=0
hnz=0
velocity=47300
totalm=3.17e-8
pnx=0           #Initial momentum 
pny=0
pnz=totalm
hnxt=hnx
hnyt=hny
hnzt=hnz
pnxt=pnx
pnyt=pny
pnzt=pnz
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
def dhdq(hi, hx, hy, hz): 
    r = np.sqrt(hx**2 + hy**2 + hz**2)
    answer =  m * hi  / r**3
    return answer

def dhdp(pi):
    answer = pi / mm
    return answer
start = time.time()
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 500:
    hnxhalf = hnx + dhdp(pnx) * dt
    hnyhalf = hny + dhdp(pny) * dt
    hnzhalf = hnz + dhdp(pnz) * dt
    pnxhalft = pnxt - dhdq(hnxt, hnxt, hnyt, hnzt) * dt
    pnyhalft = pnyt - dhdq(hnyt, hnxt, hnyt, hnzt) * dt
    pnzhalft = pnzt - dhdq(hnzt, hnxt, hnyt, hnzt) * dt
    hnx1t = hnxt + dhdp(pnxhalft) * dt * 2
    hny1t = hnyt + dhdp(pnyhalft) * dt * 2
    hnz1t = hnzt + dhdp(pnzhalft) * dt * 2
    pnx1 = pnx - dhdq(hnxhalf, hnxhalf, hnyhalf, hnzhalf) * dt * 2
    pny1 = pny - dhdq(hnyhalf, hnxhalf, hnyhalf, hnzhalf) * dt * 2
    pnz1 = pnz - dhdq(hnzhalf, hnxhalf, hnyhalf, hnzhalf) * dt * 2
    hnx1 = hnxhalf + dhdp(pnx1) * dt
    hny1 = hnyhalf + dhdp(pny1) * dt
    hnz1 = hnzhalf + dhdp(pnz1) * dt
    pnx1t = pnxhalft - dhdq(hnx1t, hnx1t, hny1t, hnz1t) * dt
    pny1t = pnyhalft - dhdq(hny1t, hnx1t, hny1t, hnz1t) * dt
    pnz1t = pnzhalft - dhdq(hnz1t, hnx1t, hny1t, hnz1t) * dt
    hnx = hnx1
    hny = hny1
    hnz = hnz1
    pnx = pnx1t
    pny = pny1t
    pnz = pnz1t
    hnxt = hnx1t
    hnyt = hny1t
    hnzt = hnz1t
    pnxt = pnx1t
    pnyt = pny1t
    pnzt = pnz1t
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