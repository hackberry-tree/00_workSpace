#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for parse_atat"""


import os
import glob
import timeit
import shutil
import unittest

import parse_atat
#import pyximport; pyximport.install(pyimport=True)
from mc import *

__date__ = "Nov 17 2014"

TEST_PATH = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
             "99_python/01_testRun/")


def test1():
    cell = CellFCC(20, mode='L10')
    print('concentration')
    print(cell.get_conc(1))
    print('local xi at 0 0 0')
    print(cell.get_tetra_xi([0, 0, 0], 0))
    print('gloval xi')
    print(cell.get_tetra_xi_glob())

lines = 'test1()'
#lines = 'cell = CellFCC(10, mode="L10")'
timeit.timeit(lines, number=1, setup="from __main__ import test1")


def main():
    unittest.main()


class Test(unittest.TestCase):
    """
    テスト
    """
    PATH = os.path.join(TEST_PATH, 'mc')

    def test_mc(self):
        """

        """
        # l = 10**5
        # print(l)
        # conc = 0
        # j = []
        # for i in range(l):
        #     j.append(int(2*round(random.uniform(conc/2., 0.5+conc/2.))-1))
        # print(j[0])
        # print(j.count(1)/float(l))
        # rand = QuadSite.random(0.5)
        # print(rand)
        # print(timeit)

        # cell = CellFCC(10)
        # print('concentration')
        # print(cell.get_conc(1))
        # print('gloval xi')
        # print(cell.get_tetra_xi_glob())

        # cell = CellFCC(10, mode='L10')
        # print('concentration')
        # print(cell.get_conc(1))
        # print('local xi at 0 0 0')
        # print(cell.get_tetra_xi([0, 0, 0], 0))
        # print('gloval xi')

        # test
        # print(cell.from_origin([0, 0, 0, 0], 1) == cell.from_origin([-1, 0, 0, 3], 2))
        # print(cell.from_origin([0, 0, 0, 3], 2) == cell.from_origin([0, 0, 0, 2], 3))
        # print(cell.from_origin([0, 0, 0, 1], 3) == cell.from_origin([0, 1, 0, 2], 0))
        # print(cell.from_origin([1, 2, 3, 1], 1) == cell.from_origin([1, 3, 4, 0], 0))
        # print(cell.from_origin([0, 0, 0, 0], 2) == cell.from_origin([0, 0, 0, 2], 0))



def clean_prev(path, files):
    """
    filesを消去
    """
    trush_list = glob.glob(os.path.join(path, files))
    for trush in trush_list:
        fname = os.path.join(path, trush)
        os.remove(fname)
        print("{0} is removed.".format(fname))

def clean_prev_dir(path, dirc):
    """
    dirctoryを消去
    """
    trush_list = glob.glob(os.path.join(path, dirc))
    for trush in trush_list:
        shutil.rmtree(trush)
        print("{0} is removed.".format(trush))


if __name__ == '__main__':
    main()
