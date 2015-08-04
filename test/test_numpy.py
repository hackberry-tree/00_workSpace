#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
numpyの記述の仕方の違いに依存したベンチマークを比較する
"""
import numpy as np
from benchmarker import Benchmarker

def for_sum(data):
    """
    for文による加算
    """
    s = 0
    for d in data:
        s += d
    return s

def builtin_sum(data):
    return sum(data)

def numpy_sum(data):
    return np.sum(data)

def alt_list_sum(data):
    """
    np.arrayをlist化してbulit inでsumを取る
    """
    return sum(list(data))

def map_for_np(data, func):
    return func(data)
vec = np.vectorize(map_for_np)

data = [1] * 2000000
data_np = np.array(data)

def bench_test_sum():
    with Benchmarker(width=20, loop=3) as bench:
        @bench("for with list")
        def _(bm):
            for_sum(data)

        @bench("builtin sum list")
        def _(bm):
            builtin_sum(data)

        @bench("np.sum numpy")
        def _(bm):
            numpy_sum(data_np)

        @bench("np.sum list")
        def _(bm):
            numpy_sum(data)

        @bench("sum numpy to list")
        def _(bm):
            alt_list_sum(data_np)

        @bench("for with numpy")
        def _(bm):
            for_sum(data_np)

        @bench("builtin sum numpy")
        def _(bm):
            builtin_sum(data_np)

loop = 200000
loop = 3
data = [1] * 200000
data_np = np.array(data)

def bench_test_for():
    """
    繰り返し処理速度を比較する
    """
    with Benchmarker(width=20, loop=3) as bench:
        @bench("for with numpy")
        def _(bm):
            l = []
            for _ in range(loop):
                l.append(numpy_sum(data_np))

        @bench("comprehension with numpy")
        def _(bm):
            l = [numpy_sum(data_np) for _ in range(loop)]

        @bench("comprehension with list")
        def _(bm):
            l = [builtin_sum(data) for _ in range(loop)]

        @bench("for with list")
        def _(bm):
            l = []
            for _ in range(loop):
                l.append(builtin_sum(data))

def bench_test_for_for():
    """
    forの入れ子
    """
    with Benchmarker(width=20, loop=3) as bench:
        @bench("plain")
        def _(bm):
            for _ in range(1000):
                numpy_sum(data_np) # 1
                numpy_sum(data_np) # 2
                numpy_sum(data_np) # 3
                numpy_sum(data_np) # 4
                numpy_sum(data_np) # 5
                numpy_sum(data_np) # 6
                numpy_sum(data_np) # 7
                numpy_sum(data_np) # 8
                numpy_sum(data_np) # 9
                numpy_sum(data_np) # 10
                numpy_sum(data_np) # 1
                numpy_sum(data_np) # 2
                numpy_sum(data_np) # 3
                numpy_sum(data_np) # 4
                numpy_sum(data_np) # 5
                numpy_sum(data_np) # 6
                numpy_sum(data_np) # 7
                numpy_sum(data_np) # 8
                numpy_sum(data_np) # 9
                numpy_sum(data_np) # 10

        @bench("for")
        def _(bm):
            for _ in range(1000):
                for _ in range(20):
                    numpy_sum(data_np) # 1

# bench_test_for()
# bench_test_sum()
bench_test_for_for()
