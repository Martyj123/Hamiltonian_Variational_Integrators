#Import useful stuffs
import math as ma
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as const

#Define constants/initial values/arrays to store stuff in
m1=1                #mass in kg of first weight
m2=1                #mass in kg of second weight
l1=1                #length in m of first stick
l2=1                #length in m of second stick
g=const.g           #Grav. acceleration in m/ss
dt=0.005             #Time step
t=0                 #Initial time


#Define the Hamiltonian of the system and the partial derivatives with respect to the different variables
def hamilt(x, y, z, u):
    a=l2**2 * m2 * z**2 + l1**2 * (m1 + m2)*u**2 - 2*m2*l1*l2*z*u*ma.cos(x-y)
    b=2*l1**2 * l2**2 * m2 *(m1 + m2*(ma.sin(x-y)**2))
    c=(m1+m2)*g*l1*ma.cos(x)
    d=m2*g*l2*ma.cos(y)
    return (a/b) - c - d

def h1(x, y, z, u):
    a=z * u * ma.sin(x - y)
    b=l1 * l2 * (m1 + m2 * ma.sin(x - y) ** 2)
    return a / b

def h2(x, y, z, u):
    a=2*m2 * l1 * l2 * z * u * ma.cos(x - y)
    b=(m1 + m2) * l1**2 * u**2
    c=m2 * (l2**2) * (z**2)
    d=(2 * (l1**2) * (l2**2) * ((m1 + m2 * (ma.sin(x - y)**2))**2))
    return (c+b-a)/d

def dhdq1(x, y, z, u):
    a=(-m1 - m2) * g * l1 * ma.sin(x)
    b=h1(x, y, z, u)
    c=h2(x, y, z, u) * ma.sin(2 * (x - y))
    return -(a-b+c)

def dhdp1(x, y, z, u):
    a=l2 * z - l1 * u * ma.cos(x - y)
    b=(l1 ** 2) * l2 * (m1 + m2 * ma.sin(x - y) ** 2)
    return a/b

def dhdq2(x, y, z, u):
    a=-m2 * g * l2 * ma.sin(y)
    b=h1(x, y, z, u)
    c=h2(x, y, z, u) * ma.sin(2 * (x - y))
    return -(a+b-c)

def dhdp2(x, y, z, u):
    a=(-m2 * l2 * z * ma.cos(x - y) + (m1 + m2) * l1 * u)
    b=m2 * l1 * (l2 ** 2) * (m1 + m2 * (ma.sin(x - y) ** 2))
    return a/b

dtvals2=[]
avgvals2=[]
#While loop to implement the algorithm
while dt < 0.05:
    q11=ma.pi/2
    p11=0
    q21=ma.pi/3
    p21=0
    qnminhalf1=q11
    pnminhalf1=p11
    qnminhalf2=q21
    pnminhalf2=p21
    evals11=[]            #store e values
    tvals=[]            #Store t values
    t=0
    while t < 100:
        q1half = qnminhalf1 + dhdp1(q11, q21, p11, p21) * dt
        p1half = pnminhalf1 - dhdq1(q11, q21, p11, p21) * dt
        q2half = qnminhalf2 + dhdp2(q11, q21, p11, p21) * dt
        p2half = pnminhalf2 - dhdq2(q11, q21, p11, p21) * dt
        qn1 = q11 + dhdp1(q1half, q2half, p1half, p2half) * dt
        pn1 = p11 - dhdq1(q1half, q2half, p1half, p2half) * dt
        qn2 = q21 + dhdp2(q1half, q2half, p1half, p2half) * dt
        pn2 = p21 - dhdq2(q1half, q2half, p1half, p2half) * dt
        #Update values for next time step
        q11 = qn1
        q21 = qn2
        p11 = pn1
        p21 = pn2
        qnminhalf1=q1half
        qnminhalf2=q2half
        pnminhalf1=p1half
        pnminhalf2=p2half
        t+=dt;
        tvals.append(t)
        E=hamilt(qn1,qn2,pn1,pn2)
        evals11.append(E)
        
    average=sum(evals11)/(100/dt)
    dtvals2.append(dt)
    avgvals2.append(average)
    print(average)
    dt=dt+0.001

#Plot the energy values
plt.plot(dtvals2, avgvals2, '-', color="black", markersize=1)
#plt.plot(tvals, p1vals, '--', color="black", markersize=1)
#plt.plot(tvals, p2vals, '-.', color="black", markersize=1)

#plt.plot(tvals, q1vals, '--', color="black", markersize=1)
#plt.plot(tvals, q2vals, '-.', color="black", markersize=1)