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
qn1=ma.pi/2       #Initial angle of the first stick
pn1=0               #Initial momentum of the first stick
qn2=ma.pi/3       #Initial angle of the second stick
pn2=0               #Initial momentum of the second stick
dt=0.1             #Time step
t=0                 #Initial time
qn1m1=qn1
pn1m1=pn1
qn2m1=qn2
pn2m1=pn2
evals21=[]            #store e values
tvals21=[]            #Store t values
posxvals1=[]        #Store position in x
posyvals1=[]        #Store position in y
posxvals2=[]        #Store position in x
posyvals2=[]        #Store position in y
q1vals=[]
q2vals=[]
p1vals=[]
p2vals=[]
error21=[]
i=0


#Define the Hamiltonian of the system and the partial derivatives with respect to the different variables
def hamilt(x, y, z, u):
    a=l2**2 * m2 * z**2 + l1**2 * (m1 + m2)*u**2 - 2*m2*l1*l2*z*u*ma.cos(x-y)
    b=2*l1**2 * l2**2 * m2 *(m1 + m2*(ma.sin(x-y)**2))
    c=(m1+m2)*g*l1*ma.cos(x)
    d=m2*g*l2*ma.cos(y)
    return (a/b) - c - d

def h1(x, y, z, u):
    a=z * u * ma.sin(x - y)
    b=l1 * l2 * (m1 + m2 * (ma.sin(x - y) ** 2))
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
    b=(l1 ** 2) * l2 * (m1 + m2 * (ma.sin(x - y) ** 2))
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

initial=hamilt(qn1, qn2, pn1, pn2)
#While loop to implement the algorithm
while t < 15:
#Verlet
    qn11 = qn1m1 + dhdp1(qn1, qn2, pn1, pn2) * dt * 2
    qn12 = qn2m1 + dhdp2(qn1, qn2, pn1, pn2) * dt* 2
    pn11 = pn1m1 - dhdq1(qn1, qn2, pn1, pn2) * dt* 2
    pn12 = pn2m1 - dhdq2(qn1, qn2, pn1, pn2) * dt* 2
    qn21 = qn1 + dhdp1(qn11, qn12, pn11, pn12) * dt* 2
    qn22 = qn2 + dhdp2(qn11, qn12, pn11, pn12) * dt* 2
    pn21 = pn1 - dhdq1(qn11, qn12, pn11, pn12) * dt* 2
    pn22 = pn2 - dhdq2(qn11, qn12, pn11, pn12) * dt* 2
    #Update values for next time step
    qn1m1 = qn11
    qn2m1 = qn12
    pn1m1 = pn11
    pn2m1 = pn12
    qn1 = qn21
    qn2 = qn22
    pn1 = pn21
    pn2 = pn22  
#Calculate/store positions for plotting
    posx1= l1 * np.sin(qn1)
    posy1= -l1 * np.cos(qn1)
    posx2= posx1 + l2 * np.sin(qn2)
    posy2= posy1 - l2 * np.cos(qn2)
    posxvals1.append(posx1)
    posyvals1.append(posy1)
    posxvals2.append(posx2)
    posyvals2.append(posy2)
    q1vals.append(qn1)
    q2vals.append(qn2)
    p1vals.append(pn1)
    p2vals.append(pn2)
#Update time step and calculate the total energy of the system at current time step
    t+=dt;
    tvals21.append(t)
    E=hamilt(qn1,qn2,pn1,pn2)
    #evals21.append((E-initial)/initial)
    evals21.append(E)
    error21.append((E-initial)/initial)
    
average=sum(evals21)/(100/dt)
print(average)
print(initial)
print(average-initial)
print(((average-initial)/initial)*100)
#Code to plot and animate the double pendulum
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def animate(i):
    thisx = [0, posxvals1[i], posxvals2[i]]
    thisy = [0, posyvals1[i], posyvals2[i]]
    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate,
                              interval=10, blit=True)


#plt.plot(tvals, q1vals, '--', color="black", markersize=1)
#plt.plot(tvals, q2vals, '-.', color="black", markersize=1)
#plt.plot(tvals, p1vals, '--', color="black", markersize=1)
#plt.plot(tvals, p2vals, '-.', color="black", markersize=1)
#Plot the energy values
#plt.subplot(1,2,1)
#plt.plot(tvals11, evals11, '-', color="black", markersize=1)
#plt.title("Scheme 2")
#plt.ylabel('Relative error')
#plt.xlabel('Time (s)')
#plt.tick_params(bottom=True, top=True, left=True, right=True)
#plt.tick_params(axis="x", direction="in", length=5, width=1)
#plt.tick_params(axis="y", direction="in", length=5, width=1)
#plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 3))
#plt.subplot(1,2,2)
#plt.plot(tvals21, error21, '-', color="black", linewidth=1, label="Relative error")
#plt.title("Scheme 1")
#plt.xlabel('Time (s)')
#plt.tick_params(bottom=True, top=True, left=True, right=True)
#plt.tick_params(axis="x", direction="in", length=5, width=1)
#plt.tick_params(axis="y", direction="in", length=5, width=1)
#plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 3))
