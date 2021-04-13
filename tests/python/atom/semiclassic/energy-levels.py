import math
import numpy as np
import matplotlib.pyplot as plt

from pyaatk.atom import SemiclassicAtom as Atom

atom = Atom(
	V    = 1.0, 
	T    = 0.01, 
	Z    = 1.0, 
	nmax = 20
)

for n in range(1,10):
	for l in range(0,n):
		enl = atom.energyLevel(n, l)
		print(n, l, enl, atom.innerRP(enl, l + 0.5), atom.outerRP(enl, l + 0.5))
