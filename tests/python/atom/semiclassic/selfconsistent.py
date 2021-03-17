import math
import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

from pyaatk.atom import SemiclassicAtom as Atom

rho = 1e-1 # g/cm^3
Avogadro = 6.022140857e+23 # N/mol
#mass = 196.966 # gold g/mol
mass = 55.845 # iron g/mol
aVol = 5.2917720859e-9**3
hartree = 13.605693009*2 # eV

V = mass/(Avogadro*rho*aVol)
T = 10/hartree # 10 eV
#Z = 79.0 # gold
Z = 26.0 # iron

atom = Atom(V=V, T=T, Z=Z, nmax=14)
r0 = atom.radius
xmax = 1.0
xmin = 1e-3
x = np.linspace(xmin, xmax, 800)**2

colors=['red', 'green', 'blue', 'orange', 'cyan', 'magenta']

grid = gs.GridSpec(1, 2)
density  = plt.subplot(grid[0,0])
potential = plt.subplot(grid[0,1])

density.plot(np.sqrt(x), atom.electronDensity(x), color="black")
potential.plot(np.sqrt(x), -r0*atom.xU(x), color="black")
print("energy levels in TF potential:")
print("n l  enl")
for n in range(1,7):
	for l in range(0,n):
		enl = atom.energyLevel(n,l)
		print(n, l, enl*hartree)

start = time.time()

Niterations = 10
print("self-consistent cycle:")
print("i  chemPot            N")
for i in range(Niterations):
	atom.update(mixing=0.75)
	print(i, atom.chemicalPotential, atom.electronStatesDiscrete())
	density.plot(np.sqrt(x), atom.electronDensity(x), color=colors[i % len(colors)])
	potential.plot(np.sqrt(x), -r0*atom.xU(x), color=colors[i % len(colors)])

elapsed = time.time() - start
print("elapsed time: ", elapsed)

print("semiclassic energy levels:")
print("n l enl")
for n in range(1,7):
	for l in range(0,n):
		enl = atom.energyLevel(n,l)
		print(n, l, enl*hartree)

plt.show()
