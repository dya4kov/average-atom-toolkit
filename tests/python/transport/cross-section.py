import math
import numpy as np
import matplotlib.pyplot as plt

from pyaatk.transport.potential import Debye, Polarization, Tabular
from pyaatk.transport import CrossSection

hbar    = 1.054571817e-34
me      = 9.10938356e-31
echarge = 1.60217662e-19
ab      = 5.29177210904e-11
ry      = 13.6056930098

ei   = 1.e-5
ef   = 1.e3
nump = 101

k = 10.0**(0.5*np.linspace(np.log10(ei/ry), np.log10(ef/ry), nump))

# potential = Debye(Z=1.0, eps=1e-4)
potential = Polarization(alpha=56.28, r0=3.0987, eps=1e-4)
var_nterms = [5, 10, 20, 30]

for nterms in var_nterms:
	print("calculating for nterms = ", nterms) 
	sigma = CrossSection(potential=potential, rmax=10000., nterms=nterms)
	cs = sigma(k)*ab*ab
	e  = ry*k*k
	plt.plot(e, cs)

plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5, 1e+3)
plt.ylim(1e-26,1e-14)
plt.show()
