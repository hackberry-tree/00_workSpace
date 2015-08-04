# -*- coding: utf-8 -*-
import numpy as np
from scipy import signal
import pylab

window = signal.gaussian(300, std=7).reshape(-1, 1)
window = np.r_[window, window][300/2:300+300/2]
pylab.plot(window)
pylab.show()

