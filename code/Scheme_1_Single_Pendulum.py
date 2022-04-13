import math as ma
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
"""
Define constants, mass, length, grav. acc., time step, initial time.
"""
m = 1; l = 10; g = 9.81; dt = 0.001; t = 0; i = 0
"""
Define initial values, position and momentum. (Position in this case is angle)
"""
qn = np.pi/4; pn = 0; qnm1 = qn; pnm1 = pn
"""
Array to store values to then place in the pandas dataframe.
"""
data = []; x = []; y = []
"""
Define the Hamiltonian of the system and its partial derivatives.
Calculate the initial Hamiltonian value of the system.
"""
def derivative_single_pend(q, p):
    h = -m*g*l*(1-ma.cos(q))+(1/2) * (p**2 / (m * l**2))
    dhdq = -m * g * l * ma.sin(q)
    dhdp = p / (m * l ** 2)
    return h, dhdq, dhdp
initial = derivative_single_pend(qn, pn)[0]
"""
Loop to iterate over the scheme for a chosen period of time.
"""
while t < 100:
    """
    Scheme one, equation (25) from the paper.
    """
    qn1 = qnm1 + derivative_single_pend(qn, pn)[2] * dt * 2
    pn1 = pnm1 - derivative_single_pend(qn, pn)[1] * dt * 2
    qn2 = qn + derivative_single_pend(qn1, pn1)[2] * dt * 2
    pn2 = pn - derivative_single_pend(qn1, pn1)[1] * dt * 2
    qn=qn2; pn=pn2; qnm1=qn1; pnm1=pn1
    """
    Update time and calculate values needed for plotting or animating.
    Store these values in the data array.
    """
    t += dt
    posx = l * ma.sin(qn)
    posy = l * ma.cos(qn)
    x.append(posx)
    y.append(posy)
    E = derivative_single_pend(qn, pn)[0]
    error = (E-initial) / initial
    data.append([t, qn, pn, posx, posy, error])
"""
Build a dataframe with the values calculated
"""
df = pd.DataFrame(data, columns = ['Time', 'Position', 'Momentum', 'X-coord', 
                                     'Y-coord', 'Error'])
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
    thisx = [0, x[i]]
    thisy = [0, y[i]]
    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate,
                              interval=0.001, blit=True)
"""
Plot the error against time for scheme 1 being tested on the single pendulum.
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