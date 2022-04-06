#Import useful stuffs
import math as ma
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as const

#Define constants/initial values/arrays to store stuff in
m1=1                #mass in kg of first weight
m2=1                #mass in kg of second weight
l1=5                #length in m of first stick
l2=7                #length in m of second stick
g=const.g           #Grav. acceleration in m/ss
q11=ma.pi/2         #Initial angle of the first stick
p11=0               #Initial momentum of the first stick
q21=ma.pi/2       #Initial angle of the second stick
p21=0               #Initial momentum of the second stick
dt=0.0001             #Time step
t=0                 #Initial time
q1half=q11
q2half=q21
p1half=p11
p2half=p21
evals=[]            #store e values
tvals=[]            #Store t values
posxvals1=[]        #Store position in x
posyvals1=[]        #Store position in y
posxvals2=[]        #Store position in x
posyvals2=[]        #Store position in y
q1vals=[]
q2vals=[]
p1vals=[]
p2vals=[]
i=0


#Define the Hamiltonian of the system and the partial derivatives with respect to the different variables
def hamilt(x, y, z, u):
    a=l2**2 * m2 * z**2 + l1**2 * (m1 + m2)*u**2 - 2*m2*l1*l2*abs(z)*abs(u)*ma.cos(x-y)
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

#While loop to implement the algorithm
while t < 100:
#Algorithm 19 from Pihajoki (Arrived at by us also)
    q1n1 = q11 + dhdp1(q1half, q2half, p1half, p2half) * dt * 2
    q2n1 = q21 + dhdp2(q1half, q2half, p1half, p2half) * dt * 2
    p1n1 = p11 - dhdq1(q1half, q2half, p1half, p2half) * dt * 2
    p2n1 = p21 - dhdq2(q1half, q2half, p1half, p2half) * dt * 2
    q1n1n = q1half + dhdp1(q1n1, q2n1, p1n1, p2n1) * dt * 2
    q2n1n = q2half + dhdp2(q1n1, q2n1, p1n1, p2n1) * dt * 2
    p1n1n = p1half - dhdq1(q1n1, q2n1, p1n1, p2n1) * dt * 2
    p2n1n = p2half - dhdq2(q1n1, q2n1, p1n1, p2n1) * dt * 2
#Update values for next time step
    q11 = q1half
    q21 = q2half
    p11 = p1half
    p21 = p2half
    q1half = q1n1n
    q2half = q2n1n
    p1half = p1n1n
    p2half = p2n1n
#Calculate/store positions for plotting
    posx1= l1 * ma.sin(q1n1n)
    posy1= -l1 * ma.cos(q1n1n)
    posx2= posx1 + l2 * ma.sin(q2n1n)
    posy2= posy1 - l2 * ma.cos(q2n1n)
    posxvals1.append(posx1)
    posyvals1.append(posy1)
    posxvals2.append(posx2)
    posyvals2.append(posy2)
    q1vals.append(q1n1n)
    q2vals.append(q2n1n)
    p1vals.append(p1n1n)
    p2vals.append(p2n1n)
#Update time step and calculate the total energy of the system at current time step
    t+=dt;
    tvals.append(t)
    E=hamilt(q1n1n,q2n1n,p1n1n,p2n1n)
    evals.append(E)

##Code to plot and animate the double pendulum
#def init():
#    line.set_data([], [])
#    time_text.set_text('')
#    return line, time_text
#
#fig = plt.figure()
#ax = fig.add_subplot(111, autoscale_on=False, xlim=(-15, 15), ylim=(-15, 15))
#ax.set_aspect('equal')
#ax.grid()
#
#line, = ax.plot([], [], 'o-', lw=2)
#time_template = 'time = %.1fs'
#time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
#
#def animate(i):
#    thisx = [0, posxvals1[i], posxvals2[i]]
#    thisy = [0, posyvals1[i], posyvals2[i]]
#    line.set_data(thisx, thisy)
#    time_text.set_text(time_template % (i*dt))
#    return line, time_text
#
#ani = animation.FuncAnimation(fig, animate,
#                              interval=1, blit=True)

##Plot the energy values
plt.plot(tvals, evals, '-', color="black", markersize=1)
#plt.plot(tvals, p1vals, '--', color="black", markersize=1)
#plt.plot(tvals, p2vals, '-.', color="black", markersize=1)

plt.plot(tvals, q1vals, '--', color="black", markersize=1)
plt.plot(tvals, q2vals, '-.', color="black", markersize=1)