import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

from pyaatk.atom import ThomasFermiAtom as Atom

grid = gs.GridSpec(1, 2)
density  = plt.subplot(grid[0,0])
potential = plt.subplot(grid[0,1])

atom = Atom(V=100.0, T=0.001, Z=100.0)
xmax = 1.0
xmin = 1e-3
x = np.linspace(xmin, xmax, 400)**2

density.plot(np.sqrt(x), 4*math.pi*atom.radius**2*x**2*atom.electronDensity(x), color="black")
potential.plot(np.sqrt(x), -atom.radius*atom.xU(x), color="black")
plt.show()