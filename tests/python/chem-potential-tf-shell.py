import numpy as np

from pyaatk.TF.eos.shell import ChemicalPotential as Mtf

M = Mtf()
M.setZ(1.0)

V = np.array([0.1, 1.0, 10.0])
T = np.array([0.1, 1.0, 10.0])

print(M(V, T))
