import math
import numpy as np
import matplotlib.pyplot as plt

from pyaatk.TF.atom import Potential as Utf
from pyaatk.transport.potential import Debye, Tabular

potential = Debye(Z=1.0)
print(potential.delta_eps(0, 100.0))
r = np.linspace(1e-6, 2.0, 50)

plt.plot(r, r*potential(r))

V = 1.0
T = 1.0
Z = 1.0
U = Utf()
U.setVTZ(V, T, Z)
r0 = (3.0*V/4.0/math.pi)**(1.0/3.0)
xmax = 1.0
xmin = 1e-3
x = np.linspace(xmin, xmax, 400)**2
r = x*r0
p = -r*U(x)
print(r0)

tfpot = Tabular(r, p)
print(tfpot.delta_eps(0, 100.0))
r = np.linspace(1e-6, 2.0, 50)
plt.plot(r, r*tfpot(r))
plt.show()
