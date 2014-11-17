#Your optional code here
#You can import some modules or create additional functions
import numpy as np
import Cython

a = [1,2,3]
b = [2,3,4]
c = [5,5,5]

d = np.array([a,b,c])

print(d.mean(axis=0))
