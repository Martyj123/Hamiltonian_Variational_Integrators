import math as ma
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
from scipy.constants import G

#Define constants/initial values/arrays to store stuff in
m=100            #mass in kg
mm=100
hnx=10           #Initial position
hny=10
hnz=10
pnx=10           #Initial momentum 
pny=10
pnz=5
hxhalf=hnx
hyhalf=hny
hzhalf=hnz
pxhalf=pnx
pyhalf=pny
pzhalf=pnz
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
    answer =  mm * hi  / r**3
    return answer

def dhdp(x):
    return x / m

#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 1000:
    hnx1 = hnx + dhdp(pxhalf) * dt
    hny1 = hny + dhdp(pyhalf) * dt
    hnz1 = hnz + dhdp(pzhalf) * dt
    pnx1 = pnx - dhdq(hxhalf, hxhalf, hyhalf, hzhalf) * dt
    pny1 = pny - dhdq(hyhalf, hxhalf, hyhalf, hzhalf) * dt
    pnz1 = pnz - dhdq(hzhalf, hxhalf, hyhalf, hzhalf) * dt
    hx2half = hxhalf + dhdp(pnx1) * dt
    hy2half = hyhalf + dhdp(pny1) * dt
    hz2half = hzhalf + dhdp(pnz1) * dt
    px2half = pxhalf - dhdq(hnx1, hnx1, hny1, hnz1) * dt
    py2half = pyhalf - dhdq(hny1, hnx1, hny1, hnz1) * dt
    pz2half = pzhalf - dhdq(hnz1, hnx1, hny1, hnz1) * dt
    hnx = hnx1
    hny = hny1
    hnz = hnz1
    pnx = pnx1
    pny = pny1
    pnz = pnz1
    pxhalf = px2half
    pyhalf = py2half
    pzhalf = pz2half
    hxhalf = hx2half
    hyhalf = hy2half
    hzhalf = hz2half
    t+=dt
    tvals.append(t)
    xpos.append(hx2half)
    ypos.append(hy2half)
    zpos.append(hz2half)
    xmom.append(px2half)
    ymom.append(py2half)
    zmom.append(pz2half)


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