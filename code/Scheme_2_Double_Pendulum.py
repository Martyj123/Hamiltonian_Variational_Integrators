import math as ma
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as const
"""
Define constants, masses, lengths, grav. acc., time step, initial time.
"""
m1=1; m2=1; l1=5; l2=7; g=const.g; dt=0.001; t=0
"""
Define initial values, position and momentum. (Position in this case is angle)
"""
q11=ma.pi/2; q21=ma.pi/2; p11=0; p21=0; q1half=q11; q2half=q21
p1half=p11; p2half=p21
"""
Array to store values to then place in the pandas dataframe.
"""
data = []; posxvals1 = []; posyvals1 = []; posxvals2 = []; posyvals2 = []
"""
Define the Hamiltonian of the system and its partial derivatives.
Calculate the initial Hamiltonian value of the system.
"""
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

def derivative_double_pend(x, y, z, u):
    a = l2**2 * m2 * z**2 + l1**2 * (m1 + m2)*u**2 - 2*m2*l1*l2*abs(z)*abs(u)*ma.cos(x-y)
    b = 2*l1**2 * l2**2 * m2 *(m1 + m2*(ma.sin(x-y)**2))
    c = (m1+m2)*g*l1*ma.cos(x)
    d = m2*g*l2*ma.cos(y)
    ham = (a/b) - c - d
    a1 = (-m1 - m2) * g * l1 * ma.sin(x)
    b1 = h1(x, y, z, u)
    c1 = h2(x, y, z, u) * ma.sin(2 * (x - y))
    dhdq1 = -(a1-b1+c1)
    a2 = l2 * z - l1 * u * ma.cos(x - y)
    b2 = (l1 ** 2) * l2 * (m1 + m2 * ma.sin(x - y) ** 2)
    dhdp1 = a2/b2
    a3 = -m2 * g * l2 * ma.sin(y)
    b3 = h1(x, y, z, u)
    c3 = h2(x, y, z, u) * ma.sin(2 * (x - y))
    dhdq2 = -(a3+b3-c3)
    a4 = (-m2 * l2 * z * ma.cos(x - y) + (m1 + m2) * l1 * u)
    b4 = m2 * l1 * (l2 ** 2) * (m1 + m2 * (ma.sin(x - y) ** 2))
    dhdp2 = a4/b4
    return ham, dhdq1, dhdq2, dhdp1, dhdp2
initial = derivative_double_pend(q11, q21, p11, p21)[0]
"""
Loop to iterate over the scheme for a chosen period of time.
"""
while t < 100:
    """
    Scheme two, equation (33) from the paper.
    """
    q1n1 = q11 + derivative_double_pend(q1half, q2half, p1half, p2half)[3] * dt
    q2n1 = q21 + derivative_double_pend(q1half, q2half, p1half, p2half)[4] * dt
    p1n1 = p11 - derivative_double_pend(q1half, q2half, p1half, p2half)[1] * dt
    p2n1 = p21 - derivative_double_pend(q1half, q2half, p1half, p2half)[2] * dt
    q1n1n = q1half + derivative_double_pend(q1n1, q2n1, p1n1, p2n1)[3] * dt
    q2n1n = q2half + derivative_double_pend(q1n1, q2n1, p1n1, p2n1)[4] * dt
    p1n1n = p1half - derivative_double_pend(q1n1, q2n1, p1n1, p2n1)[1] * dt
    p2n1n = p2half - derivative_double_pend(q1n1, q2n1, p1n1, p2n1)[2] * dt
    q11 = q1half; q21 = q2half; p11 = p1half; p21 = p2half
    q1half = q1n1n; q2half = q2n1n; p1half = p1n1n; p2half = p2n1n
    """
    Update time and calculate values needed for plotting or animating.
    Store these values in the data array.
    """
    t += dt
    posx1 = l1 * ma.sin(q1n1n)
    posy1 = -l1 * ma.cos(q1n1n)
    posx2 = posx1 + l2 * ma.sin(q2n1n)
    posy2 = posy1 - l2 * ma.cos(q2n1n)
    posxvals1.append(posx1)
    posyvals1.append(posy1)
    posxvals2.append(posx2)
    posyvals2.append(posy2)
    E = derivative_double_pend(q1n1n, q2n1n, p1n1n, p2n1n)[0]
    error = (E - initial) / initial
    data.append([t, q1n1n, q2n1n, p1n1n, p2n1n, posx1, posx2, posy1, posy2, error])
"""
Build a dataframe with the values calculated
"""
df = pd.DataFrame(data, columns = ['Time', 'Position1', 'Position2', 'Momentum1', 'Momentum2',
                                        'X-coord1', 'X-coord2', 'Y-coord1','Y-coord2', 'Error'])
"""
Animate the pendulum. The animation is super slow for some reason, changing 
the interval to be lower does nothing, can't figure it out
"""
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

fig = plt.figure(1)
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-15, 15), ylim=(-15, 15))
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
                              interval=1, blit=True)
"""
Plot the error against time for scheme 2 being tested on the double pendulum.
"""
fig2 = plt.figure(2)
plt.plot(df['Time'], df['Error'], '-', color="black", markersize=1, label="Relative error")
plt.xlabel('Time (s)')
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.ylabel('Energy')
plt.tick_params(axis="x", direction="in", length=5, width=1)
plt.tick_params(axis="y", direction="in", length=5, width=1)
plt.legend(loc='upper left', frameon=False)
fig2.show()
