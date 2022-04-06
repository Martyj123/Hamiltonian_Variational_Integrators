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
q1half=q11
p1half=p11
q21=math.pi/2       #Initial angle of the second stick
p21=0               #Initial momentum of the second stick
q2half=q21
p2half=p21
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
while t < 10:
#
    q13o2 = q1half + dhdp1(q11, q21, p11, p21) * dt
    q23o2 = q2half + dhdp2(q11, q21, p11, p21) * dt
    p13o2 = p1half/3 + dhdq1(q11, q21, p11, p21) * dt * 2/3
    p23o2 = p2half/3 + dhdq2(q11, q21, p11, p21) * dt * 2/3
    q12 = q11 + dhdp1(q13o2, q23o2, p13o2, p23o2) * dt
    q22 = q21 + dhdp2(q13o2, q23o2, p13o2, p23o2) * dt
    p12 = p11 *2/3 + dhdq1(q13o2, q23o2, p13o2, p23o2) * dt
    p22 = p21 *2/3 + dhdq2(q13o2, q23o2, p13o2, p23o2) * dt 
    q15o2 = q13o2 + dhdp1(q12, q22, p12, p22) * dt
    q25o2 = q23o2 + dhdp2(q12, q22, p12, p22) * dt
    p15o2 = p13o2 *3/2 + dhdq1(q12, q22, p12, p22) * dt * 3/2
    p25o2 = p23o2 *3/2 + dhdq2(q12, q22, p12, p22) * dt * 3/2
#Update values for next time step
    q11 = q12
    q21 = q22
    p11 = p12
    p21 = p22
    q1half = q15o2
    q2half = q25o2
    p1half = p15o2
    p2half = p25o2
#Calculate/store positions for plotting
    posx1= l1 * np.sin(q11)
    posy1= -l1 * np.cos(q11)
    posx2= posx1 + l2 * np.sin(q21)
    posy2= posy1 - l2 * np.cos(q21)
    posxvals1.append(posx1)
    posyvals1.append(posy1)
    posxvals2.append(posx2)
    posyvals2.append(posy2)
#Update time step and calculate the total energy of the system at current time step
    t+=dt;
    tvals.append(t)
    E=hamilt(q11,q21,p11,p21)
    evals.append(E)

    

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