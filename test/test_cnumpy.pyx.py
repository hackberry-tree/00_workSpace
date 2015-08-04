import numpy




def for_sum(data, s):
    for d in data:
        s += d

def builtin_sum(data):
    sum(data)

def numpy_sum(data):
    numpy.sum(data)
