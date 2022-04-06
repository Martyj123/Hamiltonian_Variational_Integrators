#Import useful stuffs
import math as math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as const

#Define constants/initial values/arrays to store stuff in
m1=1                #mass in kg of first weight
m2=1                #mass in kg of second weight
l1=2                #length in m of first stick
l2=2                #length in m of second stick
g=const.g           #Grav. acceleration in m/ss
q11=math.pi         #Initial angle of the first stick
p11=0               #Initial momentum of the first stick
q11t=q11            #Initial tilda angle of the first stick
p11t=p11            #Initial tilda momentum of the first stick
q21=math.pi/2       #Initial angle of the second stick
p21=0               #Initial momentum of the second stick
q21t=q21            #Initial tilda angle of the second stick
p21t=p21            #Initial tilda momentum of the second stick
dt=0.01             #Time step
t=0                 #Initial time
evals=[]            #store e values
tvals=[]            #Store t values
posxvals1=[]        #Store position in x
posyvals1=[]        #Store position in y
posxvals2=[]        #Store position in x
posyvals2=[]        #Store position in y
i=0


#Define the Hamiltonian of the system and the partial derivatives with respect to the different variables
def hamilt(x, y, z, u):
    a=l2**2 * m2 * u**2 + l1**2 * (m1 + m2)*u**2 - 2*m2*l1*l2*z*u*np.cos(x-y)
    b=2*l1**2 * l2**2 * m2 *(m1 + m2*np.sin(x-y)**2)
    c=(m1+m2)*g*l1*np.cos(x)
    d=m2*g*l2*np.cos(y)
    return (a/b) - c - d

def h1(x, y, z, u):
    a=z * u * np.sin(x - y)
    b=l1 * l2 * (m1 + m2 * np.sin(x - y) ** 2)
    return a / b

def h2(x, y, z, u):
    a=2*m2 * l1 * l2 * z * u * np.cos(x - y)
    b=(m1 + m2) * l1**2 * u**2
    c=m2 * (l2**2) * (z**2)
    d=(2 * (l1**2) * (l2**2) * ((m1 + m2 * (np.sin(x - y)**2))**2))
    return (c+b-a)/d

def dhdq1(x, y, z, u):
    a=(-m1 - m2) * g * l1 * np.sin(x)
    b=h1(x, y, z, u)
    c=h2(x, y, z, u) * np.sin(2 * (x - y))
    return a-b+c

def dhdp1(x, y, z, u):
    a=l2 * z - l1 * u * np.cos(x - y)
    b=(l1 ** 2) * l2 * (m1 + m2 * np.sin(x - y) ** 2)
    return a/b

def dhdq2(x, y, z, u):
    a=-m2 * g * l2 * np.sin(y)
    b=h1(x, y, z, u)
    c=h2(x, y, z, u) * np.sin(2 * (x - y))
    return a+b-c

def dhdp2(x, y, z, u):
    a=(-m2 * l2 * z * np.cos(x - y) + (m1 + m2) * l1 * u)
    b=m2 * l1 * (l2 ** 2) * (m1 + m2 * (np.sin(x - y) ** 2))
    return a/b

#While loop to implement the algorithm
while t < 100:
#Algorithm 19 from Pihajoki (Arrived at by us also)
    q1half = q11 + dhdp1(q11t, q21t, p11, p21) * dt #19a
    q2half = q21 + dhdp2(q11t, q21t, p11, p21) * dt
    p1halft = p11t + dhdq1(q11t, q21t, p11, p21) * dt #19b
    p2halft = p21t + dhdq2(q11t, q21t, p11, p21) * dt
    q1n1t = q11t + dhdp1(q1half, q2half, p1halft, p2halft) * dt * 2 #19c
    q2n1t = q21t + dhdp2(q1half, q2half, p1halft, p2halft) * dt * 2
    p1n1 = p11 + dhdq1(q1half, q2half, p1halft, p2halft) * dt * 2 #19d
    p2n1 = p21 + dhdq2(q1half, q2half, p1halft, p2halft) * dt * 2
    q1n1 = q1half + dhdp1(q1n1t, q2n1t, p1n1, p2n1) * dt #19e
    q2n1 = q2half + dhdp2(q1n1t, q2n1t, p1n1, p2n1) * dt
    p1n1t = p1halft + dhdq1(q1n1t, q2n1t, p1n1, p2n1) * dt #19f
    p2n1t = p2halft + dhdq2(q1n1t, q2n1t, p1n1, p2n1) * dt 
#Calculate/store positions for plotting
    posx1= l1 * np.sin(q1n1)
    posy1= -l1 * np.cos(q1n1)
    posx2= posx1 + l2 * np.sin(q2n1)
    posy2= posy1 - l2 * np.cos(q2n1)
    posxvals1.append(posx1)
    posyvals1.append(posy1)
    posxvals2.append(posx2)
    posyvals2.append(posy2)
#Update values for next time step
    q11 = q1n1
    q21 = q2n1
    p11 = p1n1t
    p21 = p2n1t
    q11t = q1n1
    q21t = q2n1
    p11t = p1n1t
    p21t = p2n1t
#Update time step and calculate the total energy of the system at current time step
    t+=dt;
    tvals.append(t)
    E=hamilt(q1n1,q2n1,p1n1t,p2n1t)
    evals.append(E)
##Algorithm 20 from Pihajoki (Arrived at by us also)
#    q1halft = q11t + dhdp1(q11, q21, p11t, p21t) * dt #20a
#    q2halft = q21t + dhdp2(q11, q21, p11t, p21t) * dt
#    p1half = p11 + dhdq1(q11, q21, p11t, p21t) * dt #20b 
#    p2half = p21 + dhdq2(q11, q21, p11t, p21t) * dt
#    q1n1 = q11 + dhdp1(q1halft, q2halft, p1half, p2half) * dt * 2 #20c
#    q2n1 = q21 + dhdp2(q1halft, q2halft, p1half, p2half) * dt * 2
#    p1n1t = p11t + dhdq1(q1halft, q2halft, p1half, p2half) * dt * 2 #20d
#    p2n1t = p21t + dhdq2(q1halft, q2halft, p1half, p2half) * dt * 2
#    q1n1t = q1halft + dhdp1(q1n1, q2n1, p1n1t, p2n1t) * dt #20e
#    q2n1t = q2halft + dhdp2(q1n1, q2n1, p1n1t, p2n1t) * dt
#    p1n1 = p1half + dhdq1(q1n1, q2n1, p1n1t, p2n1t) * dt #20f
#    p2n1 = p2half + dhdq2(q1n1, q2n1, p1n1t, p2n1t) * dt 
##Calculate/store positions for plotting
#    posx1= l1 * np.sin(q1n1t)
#    posy1= -l1 * np.cos(q1n1t)
#    posx2= posx1 + l2 * np.sin(q2n1t)
#    posy2= posy1 - l2 * np.cos(q2n1t)
#    posxvals1.append(posx1)
#    posyvals1.append(posy1)
#    posxvals2.append(posx2)
#    posyvals2.append(posy2)
##Update values for next time step
#    q11 = q1n1t
#    q21 = q2n1t
#    p11 = p1n1
#    p21 = p2n1
#    q11t = q1n1t
#    q21t = q2n1t
#    p11t = p1n1
#    p21t = p2n1
##Update time step and calculate the total energy of the system at current time step
#    t+=dt;
#    tvals.append(t)
#    E=hamilt(q1n1t,q2n1t,p1n1,p2n1)
#    evals.append(E)
    

#Code to plot and animate the double pendulum
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-5, 5), ylim=(-5, 5))
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

##Plot the energy values
#plt.plot(tvals[:], evals[:])