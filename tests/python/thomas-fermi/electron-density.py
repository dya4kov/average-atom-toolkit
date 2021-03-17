import math
import numpy as np
import matplotlib.pyplot as plt

from pyaatk.thomas_fermi.atom import ElectronDensity as eDens

V = 100.0
T = 0.001
Z = 100.0
rho = eDens()
rho.setVTZ(V, T, Z)
r0 = (3.0*V/4.0/math.pi)**(1.0/3.0)
xmax = 1.0
xmin = 1e-3
x = np.linspace(xmin, xmax, 400)**2
plt.plot(np.sqrt(x), 4*math.pi*r0**2*x**2*rho(x))
plt.show()