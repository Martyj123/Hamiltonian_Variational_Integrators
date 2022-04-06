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
pny=0.5
hnxt=hnx
hnyt=hny
pnxt=pnx
pnyt=pny
dt=0.1          #Time step
t=0             #Initial time
tvals2g=[]        #Store t values
xpos=[]
ypos=[]
xmom=[]
ymom=[]
evals2g=[]
error2g=[]
#Define derivatives to be used
def dhdq(hi, hx, hy): 
    r = np.sqrt(hx**2 + hy**2)
    answer =  mm * hi  / r**3
    return answer

def dhdp(pi):
    answer = pi / m
    return answer

pe=-m*mm/(np.sqrt(hnx**2+hny**2))
ke=(1/2)*(np.sqrt(pnx**2+pny**2))/m
initial = ke+pe
start = time.time()
#Loop to calculate qn2 and update qn and qn1/store the values in arrays
while t < 1000:
    hnxhalf = hnx + dhdp(pnx) * dt / 2
    hnyhalf = hny + dhdp(pny) * dt/ 2
    pnxhalft = pnxt - dhdq(hnxt, hnxt, hnyt) * dt/ 2
    pnyhalft = pnyt - dhdq(hnyt, hnxt, hnyt) * dt/ 2
    hnx1t = hnxt + dhdp(pnxhalft) * dt 
    hny1t = hnyt + dhdp(pnyhalft) * dt
    pnx1 = pnx - dhdq(hnxhalf, hnxhalf, hnyhalf) * dt 
    pny1 = pny - dhdq(hnyhalf, hnxhalf, hnyhalf) * dt 
    hnx1 = hnxhalf + dhdp(pnx1) * dt/ 2
    hny1 = hnyhalf + dhdp(pny1) * dt/ 2
    pnx1t = pnxhalft - dhdq(hnx1t, hnx1t, hny1t) * dt/ 2
    pny1t = pnyhalft - dhdq(hny1t, hnx1t, hny1t) * dt/ 2
    hnx = hnx1
    hny = hny1
    pnx = pnx1t
    pny = pny1t
    hnxt = hnx1t
    hnyt = hny1t
    pnxt = pnx1t
    pnyt = pny1t
    t+=dt
    tvals2g.append(t)
    xpos.append(hnx1)
    ypos.append(hny1)
    xmom.append(pnx)
    ymom.append(pny)
    pe=-m*mm/(np.sqrt(hnx1**2+hny1**2))
    ke=(1/2)*(np.sqrt(pnx1t**2+pny1t**2))/m
    E=ke+pe
    evals2g.append(E)
    error2g.append((E-initial)/initial)
        
average=sum(evals2g)/(t/dt)
print(average)
print(initial)
print(average-initial)
print(((average-initial)/initial)*100)

end = time.time()
#plt.plot(tvals,ymom,'.', color="black", markersize=1)
plt.plot(tvals2g, evals2g, '.', color="black", markersize=1)
plt.xlabel('x')
plt.ylabel('y')