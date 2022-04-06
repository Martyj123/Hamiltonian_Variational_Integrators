def h1(x, y, z, u):
    return ((z * u * np.sin(x - y)) / (l1 * l2 * (m1 + m2 * np.sin(x - y) * np.sin(x - y))))

def h2(x, y, z, u):
    return ((m2 * l2**2 * z**2 + (m1 + m2) * l1**2 * u**2 - 2 * m2 * l1 * l2 * z * u * np.cos(x - y)) / (2 * l1**2 * l2**2 * (m1 + m2 * np.sin(x - y) * np.sin(x - y))**2))

def dhdq1(x, y, z, u):
    return (-(m1 + m2) * g * l1 * np.sin(x) - h1(x, y, z, u) + h2(x, y, z, u) * np.sin(2 * (x - y)))

def dhdp1(x, y, z, u):
    return ((l2 * z - l1 * u * np.cos(x - y)) / (l1 **2 * l2 * (m1 + m2 * np.sin(x - y) * np.sin(x - y))))

def dhdq2(x, y, z, u):
    return (-m2 * g * l2 * np.sin(y) + h1(x, y, z, u) - h2(x, y, z, u) * np.sin(2 * (x - y)))

def dhdp2(x, y, z, u):
    return ((-m2 * l2 * z * np.cos(x - y) + (m1 + m2) * l1 * u) / (m2 * l1 * l2 **2 * (m1 + m2 * np.sin(x - y) * np.sin(x - y))))