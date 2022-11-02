import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import RK45

# Data :
graph = True
L = 3  # Length [m]
d = 0.03  # Diameter [m]
p = 2  # Pressure [bar]
deltaH = 20  # Enthalpy [J/kmol]
area = np.pi * (d / 2) ** 2
MA = 50  # Molar mass [kg/kmol]
cp = 1.2  # Heat Capacity [J/kg.K]
Tf = 650  # Furnace temperature [K]
U = 140  # Heat transfer coefficient[W/m2*K]
mdot = MA
n = 100
A = np.pi * d / area

# Explicit Euler method

h = (3 - 0) / n
Z = np.linspace(0, 3, n + 1)
T = np.zeros(n + 1)
X = np.zeros(n + 1)
X[0] = 0
T[0] = 500
X1 = np.zeros(n + 1)
T1 = np.zeros(n + 1)
X1[0] = 0
T1[0] = 500

for i in range(n):

    r = 4 * (10 ** 5) * (np.exp(-25000 / (8.314 * T[i]))) * 2 * (1 - X[i])   # [kmol/(m^2 * s)]
    X[i + 1] = X[i] + h * r * area
    print(X[i])
    T[i + 1] = T[i] + h * r * area * (-deltaH) * 1000 / (mdot * cp)

    r1 = 4 * (10 ** 5) * (np.exp(-25000 / (8.314 * T1[i]))) * 2 * (1 - X1[i])
    X1[i + 1] = X1[i] + h * r * area
    T1[i + 1] = T1[i] + h * (r * (-deltaH) * area / (mdot * cp) + U * 3600 * A * area * (Tf - T1[i]) / (cp * mdot))

#
# Graphical part
#
if graph is True:
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(5, 10), layout='tight')
    fig.suptitle('Plug Flow Reactor Characteristics')

    ax1.plot(Z, X, 'b')
    ax1.plot(Z, X1, 'r')
    ax1.set_title('Conversion graph')
    ax1.set(xlabel='Length [m]', ylabel='Conversion [%]')
    ax2.plot(Z, T, 'b')
    ax2.plot(Z, T1, 'r')
    ax2.set_title('Temperature graph')
    ax2.set(xlabel='Length [m]', ylabel='Temperature [K]')
    plt.show()
