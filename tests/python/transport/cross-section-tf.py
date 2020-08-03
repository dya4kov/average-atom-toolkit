import math
import numpy as np
import matplotlib.pyplot as plt

from pyaatk.TF.atom import Potential as Utf
from pyaatk.transport.potential import Tabular
from pyaatk.transport import CrossSection

hbar    = 1.054571817e-34
me      = 9.10938356e-31
echarge = 1.60217662e-19
ab      = 5.29177210904e-11
ry      = 13.6056930098

ei   = 1.e-5
ef   = 1.e3
nump = 101

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
tfpot = Tabular(r, p, eps=1e-4)

k = 10.0**(0.5*np.linspace(np.log10(ei/ry), np.log10(ef/ry), nump))

sigma = CrossSection(potential=tfpot, rmax=10000., nterms=30)
cs = sigma(k)*ab*ab
e  = ry*k*k

plt.plot(e, cs)

plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5, 1e+3)
plt.ylim(1e-26,1e-14)
plt.show()
