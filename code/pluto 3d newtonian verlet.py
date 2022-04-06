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
m=100           #mass in kg
mm=100
hnx=-10         #Distance at aphelion in km
hny=0
velocity=47360/c           #min orbital velocity in km/s
pnx=10          #Initial momentum
pny=5
hnxt=hnx
hnyt=hny
pnxt=pnx
pnyt=pny
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
def dhdq(hi, hx, hy): 
    r = np.sqrt(hx**2 + hy**2)
    answer =  m * hi  / r**3
    return answer

def dhdp(pi):
    answer = pi / mm
    return answer
start = time.time()
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 50:
    hnxhalf = hnx + dhdp(pnx) * dt
    hnyhalf = hny + dhdp(pny) * dt
    pnxhalft =pnxt - dhdq(hnxt, hnxt, hnyt) * dt
    pnyhalft = pnyt - dhdq(hnyt, hnxt, hnyt) * dt
    hnx1t = hnxt + dhdp(pnxhalft) * dt * 2
    hny1t = hnyt + dhdp(pnyhalft) * dt * 2
    pnx1 = pnx - dhdq(hnxhalf, hnxhalf, hnyhalf) * dt * 2
    pny1 = pny - dhdq(hnyhalf, hnxhalf, hnyhalf) * dt * 2
    hnx1 = hnxhalf + dhdp(pnx1) * dt
    hny1 = hnyhalf + dhdp(pny1) * dt
    pnx1t = pnxhalft - dhdq(hnx1t, hnx1t, hny1t) * dt
    pny1t = pnyhalft - dhdq(hny1t, hnx1t, hny1t) * dt
    hnx = hnx1
    hny = hny1
    pnx = pnx1t
    pny = pny1t
    hnxt = hnx1t
    hnyt = hny1t
    pnxt = pnx1t
    pnyt = pny1t
    t+=dt
    tvals.append(t)
    xpos.append(hnx1)
    ypos.append(hny1)
    xmom.append(pnx)
    ymom.append(pny)

end = time.time()
print(end - start)
plt.plot(xpos, ypos)
plt.xlabel('x')
plt.ylabel('y')