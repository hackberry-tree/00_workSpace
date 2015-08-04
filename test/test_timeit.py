import timeit

def timeit_list(n, iter):
    list_setup = """
import test_numpy
data = [1] * {}
s = 0
""".format(n)
    print('list(python):', end=' ')
    print(timeit.timeit(
        "test_numpy.for_sum(data,s)", list_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_numpy.builtin_sum(data)", list_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_numpy.numpy_sum(data)", list_setup, number=iter))

def timeit_list_c(n, iter):
    list_setup = """
import pyximport; pyximport.install(pyimport=True)
import test_cnumpy
data = [1] * {}
s = 0
""".format(n)
    print('list(cython):', end=' ')
    print(timeit.timeit(
        "test_cnumpy.for_sum(data,s)", list_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_cnumpy.builtin_sum(data)", list_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_cnumpy.numpy_sum(data)", list_setup, number=iter))

def timeit_array(n, iter):
    array_setup = """
import numpy
import array
import test_numpy
data = array.array('L', [1] * {})
s = 0
""".format(n)
    print('array(python):', end=' ')
    print(timeit.timeit(
        "test_numpy.for_sum(data,s)", array_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_numpy.builtin_sum(data)", array_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_numpy.numpy_sum(data)", array_setup, number=iter))

def timeit_array_c(n, iter):
    array_setup = """
import numpy
import array
import pyximport; pyximport.install(pyimport=True)
import test_cnumpy
data = array.array('L', [1] * {})
s = 0
""".format(n)
    print('array(cython):', end=' ')
    print(timeit.timeit(
        "test_cnumpy.for_sum(data,s)", array_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_cnumpy.builtin_sum(data)", array_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_cnumpy.numpy_sum(data)", array_setup, number=iter))

def timeit_numpy(n, iter):
    numpy_setup = ""
import numpy
import test_numpy
data = numpy.array([1] * {})
s = 0
""".format(n)
    print('numpy.array(python):', end=' ')
    print(timeit.timeit(
        "test_numpy.for_sum(data,s)", numpy_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_numpy.builtin_sum(data)", numpy_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_numpy.numpy_sum(data)", numpy_setup, number=iter))

def timeit_numpy_c(n, iter):
    numpy_setup = """
import numpy
import pyximport; pyximport.install(pyimport=True)
import test_cnumpy
data = numpy.array([1] * {})
s = 0
""".format(n)
    print('numpy.array(cython):', end=' ')
    print(timeit.timeit(
        "test_cnumpy.for_sum(data,s)", numpy_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_cnumpy.builtin_sum(data)", numpy_setup, number=iter), end=' ')
    print(timeit.timeit(
        "test_cnumpy.numpy_sum(data)", numpy_setup, number=iter))


if __name__ == '__main__':
    #timeit.timeit("test_numpy.timeit_list(1,1)", 'import test_numpy', number=1)
    timeit_list(100000, 1000)
    timeit_list_c(100000, 1000)
    timeit_array(100000, 1000)
    timeit_array_c(100000, 1000)
    timeit_numpy(100000, 1000)
    timeit_numpy_c(100000, 1000)
