import math as math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#Define constants/initial values/arrays to store stuff in
m1=1             #mass in kg
m2=1             #mass in kg
l1=1            #length in m
l2=1            #length in m
g=9.81         #Grav. acceleration in m/ss
q11=math.pi/2       #Initial angle of the pendulum
p11=1            #Initial momentum of the pendulum
q12=math.pi/2     #Angle at first time step
p12=1           #Momentum at first time step
q21=math.pi/3        #Initial angle of the pendulum
p21=1            #Initial momentum of the pendulum
q22=math.pi/3     #Angle at first time step
p22=1           #Momentum at first time step
dt=0.01         #Time step
t=0             #Initial time
qvals=[]        #Store q values
pvals=[]        #Store p values
tvals=[]        #Store t values
posxvals1=[]     #Store position in x
posyvals1=[]     #Store position in y
posxvals2=[]     #Store position in x
posyvals2=[]     #Store position in y
n=0
def h1(x, y, z, u):
    return ((z * u * math.sin(x - y)) / (l1 * l2 * (m1 + m2 * math.sin(x - y) ** 2)))

def h2(x, y, z, u):
    a=m2 * l1 * l2 * z * u * math.cos(x - y)
    b=(m1 + m2) * l1**2 * u**2
    c=m2 * (l2**2) * (z**2)
    d=(2 * (l1**2) * (l2**2) * ((m1 + m2 * (math.sin(x - y)**2))**2))
    return (c+b-a)/d

def dhdq1(x, y, z, u):
    return ((-m1 - m2) * g * l1 * math.sin(x) - h1(x, y, z, u) + h2(x, y, z, u) * math.sin(2 * (x - y)))

def dhdp1(x, y, z, u):
    return ((l2 * z - l1 * u * math.cos(x - y)) / ((l1 ** 2) * l2 * (m1 + m2 * math.sin(x - y) ** 2)))

def dhdq2(x, y, z, u):
    return (-m2 * g * l2 * math.sin(y) + h1(x, y, z, u) - h2(x, y, z, u) * math.sin(2 * (x - y)))

def dhdp2(x, y, z, u):
    return ((-m2 * l2 * z * math.cos(x - y) + (m1 + m2) * l1 * u) / (m2 * l1 * (l2 ** 2) * (m1 + m2 * (math.sin(x - y) ** 2))))

while t < 10:
    qn1 = 2 * dhdp1(q12, q22, p12, p22) * dt + q11
    print("qn1 is ", qn1)
    pn1 = p11 - 2 * -dhdq1(q12, q22, p12, p22) * dt
    print("pn1 is ", pn1)
    qn2 = 2 * dhdp2(q12, q22, p12, p22) * dt + q21
    print("qn2 is ", qn2)
    pn2 = p21 - 2 * -dhdq2(q12, q22, p12, p22) * dt
    print("pn2 is ", pn2)
    print(n)
    n+=1
    q11=q12
    p11=p12
    q21=q22
    p21=p22
    q12=qn1
    p12=pn1
    q22=qn2
    p22=pn2
    t+=dt;
    tvals.append(t)
    posx1= l1 * math.sin(qn1)
    posy1= -l1 * math.cos(qn1)
    posx2= posx1 + l2 * math.sin(qn2)
    posy2= posy1 - l2 * math.cos(qn2)
    posxvals1.append(posx1)
    posyvals1.append(posy1)
    posxvals2.append(posx2)
    posyvals2.append(posy2)
    
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