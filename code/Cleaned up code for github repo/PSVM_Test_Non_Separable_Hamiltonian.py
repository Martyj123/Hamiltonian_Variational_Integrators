import math as ma
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
"""
Define constants: time step, initial time.
"""
dt = 0.001; t = 0; i = 0
"""
Define initial values, position and momentum. (Position in this case is angle)
"""
qn = -3; pn = 0
"""
Array to store values to then place in the pandas dataframe.
"""
data = []
"""
Define the Hamiltonian of the system and its partial derivatives.
Calculate the initial Hamiltonian value of the system.
"""
def derivative_test_hamilt(q, p):
    h = (1/2)*(q**2*p**2 +q**2+(p**2)+1)
    dhdq = q*p**2 + q
    dhdp = p*q**2 + p
    return h, dhdq, dhdp
initial = derivative_test_hamilt(qn, pn)[0]
"""
Loop to iterate over the scheme for a chosen period of time.
"""
while t < 100:
    """
    PSVM scheme, equation (29) from the paper.
    """
    qnh = qn + derivative_test_hamilt(qn, pn)[2] * dt / 2
    pnh = pn - derivative_test_hamilt(qn, pn)[1] * dt / 2
    qn1 = qn + derivative_test_hamilt(qnh, pnh)[2] * dt
    pn1 = pn - derivative_test_hamilt(qnh, pnh)[1] * dt
    qn2 = qnh + derivative_test_hamilt(qn1, pn1)[2] * dt / 2
    pn2 = pnh - derivative_test_hamilt(qn1, pn1)[1] * dt / 2
    qn=qn2; pn=pn2
    """
    Update time and calculate values needed for plotting or animating.
    Store these values in the data array.
    """
    t += dt
    E = derivative_test_hamilt(qn, pn)[0]
    error = (E-initial) / initial
    data.append([t, qn, pn, error])
"""
Build a dataframe with the values calculated
"""
df = pd.DataFrame(data, columns = ['Time', 'Position', 'Momentum', 'Error'])
"""
Plot the error against time for PSVM scheme being tested on the non-separable
Hamiltonian.
"""
fig2 = plt.figure()
plt.plot(df['Time'], df['Error'], '-', color="black", markersize=1, label="Relative error")
plt.xlabel('Time (s)')
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.ylabel('Energy')
plt.tick_params(axis="x", direction="in", length=5, width=1)
plt.tick_params(axis="y", direction="in", length=5, width=1)
plt.legend(loc='upper left', frameon=False)
fig2.show()