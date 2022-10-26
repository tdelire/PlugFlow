import numpy as np
from matplotlib import pyplot as plt

# Data :
L = 3  # Length [m]
d = 0.2  # Diameter [m]
p = 2  # Pressure [bar]
T0 = 500  # Temperature [K]
k = 4 * (10 ** 5) * np.exp(-25000 / (8.314 * T0))
r = k * p  # Reaction rate [mol/s]
deltaH = 20000  # Enthalpy [J/kmol]
area = np.pi * (d / 2) ** 2
MA = 50  # Molar mass [kg/kmol]
cp = 1200  # Heat Capacity [J/kg.K]
Tf = 650  # Furnace temperature [K]
U = 140  # Heat transfer coefficient[W/m2*K]
mdot = 1000 * MA / 3600
# Explicit Euler method

Zstart = 0
Zend = 3
Tstart = 500
Xstart = 0

for n in [25]:
    h = (Zend - Zstart) / n
    Z = np.linspace(Zstart, Zend, n + 1)
    T = np.zeros(n + 1)
    X = np.zeros(n + 1)
    X[0] = Xstart
    T[0] = Tstart
    for i in range(n):
        k = 4 * (10 ** 5) * np.exp(-25000 / (8.314 * T[i]))
        r = k * p / 50    # [kmol/(m^2 * s)]
        X[i + 1] = X[i] + h * Z[i] * r * area * 1000 / 3600
        T[i + 1] = T[i] + h * Z[i] * r * area * (-deltaH) / (mdot * cp)

#
# Graphical part
#

fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(5, 10), layout='tight')
fig.suptitle('Plug Flow Reactor Characteristics')

ax1.plot(Z, X)
ax1.set_title('Conversion graph')
ax1.set(xlabel='Length [m]', ylabel='Conversion [/]')
ax2.plot(Z, T)
ax2.set_title('Temperature graph')
ax2.set(xlabel='Length [m]', ylabel='Temperature [K]')
plt.show()
