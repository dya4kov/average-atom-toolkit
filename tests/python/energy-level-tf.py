import math
import numpy as np
import matplotlib.pyplot as plt

from pyaatk.thomas_fermi.atom import EnergyLevel as eLevel
from pyaatk.thomas_fermi.atom import RotatePoints as RP

V = 1.0
T = 1.0
Z = 100.0
e = eLevel()
rp = RP()
e.setVTZ(V, T, Z)
rp.setVTZ(V, T, Z)

for n in range(1,7):
	for l in range(0,n):
		enl = e(n,l)
		print(n, l, enl, rp.inner(enl, l + 0.5), rp.outer(enl, l + 0.5))
