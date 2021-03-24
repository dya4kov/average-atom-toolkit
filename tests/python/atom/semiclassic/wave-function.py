import math
import numpy as np
import matplotlib.pyplot as plt

from mendeleev import element

from pyaatk.atom import SemiclassicAtom as Atom

# rho = 1e-1 # g/cm^3
Avogadro = 6.022140857e+23 # N/mol
aVol = 5.2917720859e-9**3
hartree = 13.605693009*2 # eV

elem     = element('Al')
rho      = elem.density   #2.7 g/cm^3
mass     = elem.atomic_weight # 27 g/mol
# mass = 196.966 # gold g/mol
# mass = 55.845 # iron g/mol
Z        = elem.atomic_number #13.0
# Z = 79.0 # gold
# Z = 26.0 # iron

V = mass/(Avogadro*rho*aVol)
T = 77/hartree # 10 eV

atom = Atom(V=V, T=T, Z=Z, nmax=14)
r0 = (3.0*V/4.0/math.pi)**(1.0/3.0)
xmax = 1.0
xmin = 1e-3
x = np.linspace(xmin, xmax, 1500)**2
atom.update(0.75);

for n in range(1, 3):
	for l in range(0,n):
		enl = atom.energyLevel(n, l)
		Rnl = atom.waveFunction(x, enl, l + 0.5)
		plt.plot(np.sqrt(x), Rnl)
		print(n, l, enl, atom.innerRP(enl, l + 0.5), atom.outerRP(enl, l + 0.5))

# plt.xlim(0, 0.25)
plt.show()
