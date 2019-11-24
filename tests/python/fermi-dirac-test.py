import numpy as np
import matplotlib.pyplot as plt

from pynumtk.specfunc.FermiDirac import ThreeHalf as FD3half
from pynumtk.specfunc.FermiDirac import Half as FDhalf
from pynumtk.specfunc.FermiDirac import MHalf as FDmhalf

x = np.linspace(-1.0, 4.0, 20)

f1 = FD3half()
f2 = FDhalf()
f3 = FDmhalf()

plt.plot(x, f1(x))
plt.plot(x, f2(x))
plt.plot(x, f3(x))
plt.show()