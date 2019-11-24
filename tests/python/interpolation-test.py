import math
import numpy as np
import matplotlib.pyplot as plt

from pynumtk.interpolation import Polynomial as Interp


x = np.linspace(1.0, 20.0, 20)
y = np.sin(x)

yi = Interp(x, y, 3)
xi = np.linspace(1.0, 20.0, 100)

plt.plot(x, y)
plt.plot(xi, yi(xi))

plt.show()