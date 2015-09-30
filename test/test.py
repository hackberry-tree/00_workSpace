# -*- coding: utf-8 -*-
import numpy as np
from scipy import signal
import pylab

spins = np.array([-1, 0, 1, 2, 4])
print(spins.shape)

pos = np.array([np.array([1, 2, 3]), np.array([2, 3])])
print(pos.prod(0))
print(pos)
print(pos.shape)
xis = spins[pos]
print(xis)
