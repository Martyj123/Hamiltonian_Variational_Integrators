import math as ma
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
from scipy.constants import G
import time
#Define constants/initial values/arrays to store stuff in
m=10             #mass in kg
mm=10
hnx=10           #Initial position
hny=10
pnx=2           #Initial momentum 
pny=0
hnxminhalf=hnx
hnyminhalf=hny
pnxminhalf=pnx
pnyminhalf=pny
dt=0.001          #Time step
t=0             #Initial time
tvals=[]        #Store t values
xpos=[]
ypos=[]
xmom=[]
ymom=[]
evalsn2=[]
errorn2=[]
#Define derivatives to be used
def dhdq(hi, hx, hy): 
    r = np.sqrt(hx**2 + hy**2)
    answer =  mm * hi  / r**3
    return answer

def dhdp(pi):
    answer = pi / m
    return answer
pe=-m*mm/(ma.sqrt(hnx**2+hny**2)*2)
ke=(1/2)*(np.sqrt(pnx**2+pny**2))/m
initial = pe
start = time.time()
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 200:
    hnxhalf = hnxminhalf + dhdp(pnx) * dt
    hnyhalf = hnyminhalf + dhdp(pny) * dt
    pnxhalf = pnxminhalf - dhdq(hnx, hnx, hny) * dt
    pnyhalf = pnyminhalf - dhdq(hny, hnx, hny) * dt
    hnx1 = hnx + dhdp(pnxhalf) * dt 
    hny1 = hny + dhdp(pnyhalf) * dt
    pnx1 = pnx - dhdq(hnxhalf, hnxhalf, hnyhalf) * dt 
    pny1 = pny - dhdq(hnyhalf, hnxhalf, hnyhalf) * dt 
    hnx = hnx1
    hny = hny1
    pnx = pnx1
    pny = pny1
    hnxminhalf = hnxhalf
    hnyminhalf = hnyhalf
    pnxminhalf = pnxhalf
    pnyminhalf = pnyhalf
    t+=dt
    tvals.append(t)
    xpos.append(hnx1)
    ypos.append(hny1)
    xmom.append(pnx)
    ymom.append(pny)
    pe=-m*mm/(np.sqrt(hnx1**2+hny1**2)*2)
    ke=(1/2)*(np.sqrt(pnx1**2+pny1**2))/m
    E=pe
    evalsn2.append(E)
    errorn2.append((E-initial)/initial)
        
average=sum(evalsn2)/(200/dt)
print(average)
print(initial)
print(average-initial)
print((average-initial)/initial*100)

end = time.time()
#plt.plot(tvals,ymom,'.', color="black", markersize=1)
plt.plot(tvals, evalsn2, '.', color="black", markersize=1)
plt.xlabel('x')
plt.ylabel('y')