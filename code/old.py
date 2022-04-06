from numpy import array,arange,add
import numpy as np
from pylab import plot,xlabel,ylabel,show,figure

#Defining constants
G=1 #m^3kg^-1s^-2
M=1 #M
m1=0.000954 #M
m2=2.86e-4 #M
AU=1.495978707e11 #m
a1=5.204 #AU
a2=9.583 #AU
#Defining the function governing forces between the objects
#1 representing first object, 2 the second
#10 distance between obj. 1 and the sun, 20 obj. 2 to the sun
#12 distance between the two objects
def func(r,t):
    posx1=r[0]
    posy1=r[1]
    velx1=r[2]
    vely1=r[3]
    posx2=r[4]
    posy2=r[5]
    velx2=r[6]
    vely2=r[7]
    dposx1dt=velx1
    dposy1dt=vely1
    dposx2dt=velx2
    dposy2dt=vely2
    rval10=np.sqrt(posx1**2+posy1**2)
    rval20=np.sqrt(posx2**2+posy2**2)
    rval12=np.sqrt((posx1-posx2)**2+(posy1-posy2)**2)
    dvelx1dt=(-(G*M*(posx1)/rval10**3)+(G*m2*(posx1-posx2)/rval12**3))
    dvely1dt=(-(G*M*(posy1)/rval10**3)+(G*m2*(posy1-posy2)/rval12**3))
    dvelx2dt=(-(G*M*(posx2)/rval20**3)-(G*m1*(posx1-posx2)/rval12**3))
    dvely2dt=(-(G*M*(posy2)/rval20**3)-(G*m1*(posy1-posy2)/rval12**3))
    return array([dposx1dt,dposy1dt,dvelx1dt,dvely1dt,dposx2dt,dposy2dt,dvelx2dt,dvely2dt],float)
#Initialising the tpoints used in RK4
a = 0
b = 160
N = 1000
h = (b-a)/N
tpoints = arange(a,b,h)
r=array([a1,0,0,0.4,a2,0,0,0.3],float) #Positions in m, velocities in m/s
posx1_pnts = []
posy1_pnts = []
posx2_pnts = []
posy2_pnts = []
velx1_pnts = []
vely1_pnts = []
velx2_pnts = []
vely2_pnts = []
k1=[]*4
k2=[]*4
k3=[]*4
k4=[]*4
for t in tpoints:
    k1 = h*func(r,t)
    k2 = h*func(r+0.5*k1,t+0.5*h)
    k3 = h*func(r+0.5*k2,t+0.5*h)
    k4 = h*func(r+k3,t+h)
    r = add(r, (k1+2*k2+2*k3+k4)/6)
    posx1_pnts.append(r[0])
    posy1_pnts.append(r[1])
    velx1_pnts.append(r[2])
    vely1_pnts.append(r[3])
    posx2_pnts.append(r[4])
    posy2_pnts.append(r[5])
    velx2_pnts.append(r[6])
    vely2_pnts.append(r[7])
#1 is in blue 2 is in orange
figure(num=None, figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')
plot(posx1_pnts,posy1_pnts,linestyle="",marker=".",markersize=1)
plot(posx2_pnts,posy2_pnts,linestyle="",marker=".",markersize=1)
plot(0,0,linestyle="",marker=".",markersize=5)
xlabel('x [AU]')
ylabel('y [AU]')
show()
