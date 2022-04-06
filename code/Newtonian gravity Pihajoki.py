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
dt=0.01          #Time step
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
    answer =  mm * hi  / r**3
    return answer

def dhdp(pi):
    answer = pi / m
    return answer
start = time.time()
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 300:
    hnxhalf = hnx + dhdp(pnx) * dt / 2
    hnyhalf = hny + dhdp(pny) * dt/ 2
    hnzhalf = hnz + dhdp(pnz) * dt/ 2
    pnxhalft = pnxt - dhdq(hnxt, hnxt, hnyt, hnzt) * dt/ 2
    pnyhalft = pnyt - dhdq(hnyt, hnxt, hnyt, hnzt) * dt/ 2
    pnzhalft = pnzt - dhdq(hnzt, hnxt, hnyt, hnzt) * dt/ 2
    hnx1t = hnxt + dhdp(pnxhalft) * dt 
    hny1t = hnyt + dhdp(pnyhalft) * dt
    hnz1t = hnzt + dhdp(pnzhalft) * dt 
    pnx1 = pnx - dhdq(hnxhalf, hnxhalf, hnyhalf, hnzhalf) * dt 
    pny1 = pny - dhdq(hnyhalf, hnxhalf, hnyhalf, hnzhalf) * dt 
    pnz1 = pnz - dhdq(hnzhalf, hnxhalf, hnyhalf, hnzhalf) * dt 
    hnx1 = hnxhalf + dhdp(pnx1) * dt/ 2
    hny1 = hnyhalf + dhdp(pny1) * dt/ 2
    hnz1 = hnzhalf + dhdp(pnz1) * dt/ 2
    pnx1t = pnxhalft - dhdq(hnx1t, hnx1t, hny1t, hnz1t) * dt/ 2
    pny1t = pnyhalft - dhdq(hny1t, hnx1t, hny1t, hnz1t) * dt/ 2
    pnz1t = pnzhalft - dhdq(hnz1t, hnx1t, hny1t, hnz1t) * dt/ 2
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