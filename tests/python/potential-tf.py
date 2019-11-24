import math
import numpy as np
import matplotlib.pyplot as plt

from pyaatk.TF.atom import Potential as Utf

V = 1.0
T = 1.0
Z = 1.0
U = Utf()
U.setVTZ(V, T, Z)
r0 = (3.0*V/4.0/math.pi)**(1.0/3.0)
xmax = 1.0
xmin = 1e-3
x = np.linspace(xmin, xmax, 400)**2
plt.plot(np.sqrt(x), -x*r0*U(x))
plt.show()