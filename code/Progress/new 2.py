import math as ma
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#Define constants/initial values/arrays to store stuff in
m1=1             #mass in kg
m2=1             #mass in kg
l1=10            #length in m
l2=10            #length in m
g=-9.81         #Grav. acceleration in m/ss
q11=np.pi/4        #Initial angle of the pendulum
p11=0            #Initial momentum of the pendulum
q12=np.pi/4     #Angle at first time step
p12=0           #Momentum at first time step
q21=np.pi/4        #Initial angle of the pendulum
p21=0            #Initial momentum of the pendulum
q22=np.pi/4     #Angle at first time step
p22=0           #Momentum at first time step
dt=0.01         #Time step
t=0             #Initial time
qvals=[]        #Store q values
pvals=[]        #Store p values
tvals=[]        #Store t values
posxvals=[]     #Store position in x
posyvals=[]     #Store position in y

h1= (p12 * p22 * np.sin(q12 - q22)) / (2 * l1 * l2 * (m1 + m2 * np.sin(q12 - q22) * np.sin(q12 - q22))**2)

h2= (m2 * l2**2 * p12**2 + (m1 + m2) * l1**2 * p22**2 - 2 * m2 * l1 * l2 * p11 * p22 * np.cos(q12 - q22))/ (2 * l1**2 * l2**2 * (m1 + m2 * np.sin(q12 - q22) * np.sin(q12 - q22))**2)

def dhdq1(x, y, z, u):
    return -(m1 + m2)*g*l1*np.sin(x) - h1 + h2 * np.sin(2*(x-y))

def dhdp1(x, y, z, u):
    return (l2 * p12 - l1 * p22 * np.cos(q12 - q22)) / (l1 ** 2 * l2 * (m1 + m2 * np.sin(q12 - q22) * np.sin(q12 - q22)))

def dhdq2(x, y, z, u):
    return - m2 * g * l2 * np.sin(q22) + h1 - h2 * np.sin(2 * (q12 - q22))

def dhdp2(x, y, z, u):
    return (-m2 * l2 * p12 * np.cos(q12 - q22) + (m1 + m2) * l1 * p22) / (m2 * l1 * l2 ** 2 * (m1 + m2 * np.sin(q12 - q22) * np.sin(q12 - q22)))

while t < 20:
    qn1 = 2 * dhdp(pn1) * dt + qn
    pn1 = pn - 2 * dhdq(qn1) * dt
    qn2 = 2 * dhdp(pn1) * dt + qn
    pn2 = pn - 2 * dhdq(qn1) * dt
    qvals.append(qn2)
    pvals.append(pn2)
    qn=qn1
    pn=pn1
    qn1=qn2
    pn1=pn2
    t+=dt;
    tvals.append(t)
    posx=l * ma.sin(qn2)
    posy=l * ma.cos(qn2)
    posxvals.append(posx)
    posyvals.append(posy)