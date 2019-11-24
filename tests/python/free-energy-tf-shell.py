import numpy as np

from pyaatk.TF.eos.shell import FreeEnergy as Ftf

F = Ftf()
F.setZ(1.0)
F.setNmax(25)

V = np.array([0.1, 1.0, 10.0])
T = np.array([0.1, 1.0, 10.0])

print(F(V, T))
print(F.DV(V, T))
print(F.DT(V, T))
print(F.D2V(V, T))
print(F.DVT(V, T))
print(F.D2T(V, T))
