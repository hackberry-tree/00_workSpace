#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for parse_atat"""


import os
import glob
import shutil
import unittest
from mc_fct import *
# import numpy as np
# cimport numpy as np

__date__ = "Nov 17 2014"

TEST_PATH = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
             "99_python/01_testRun/")


def main():
    unittest.main()


class Test(unittest.TestCase):
    """
    テスト
    """
    PATH = os.path.join(TEST_PATH, 'mc')

    def _test_tetrahedron(self):
        xis = Tetrahedron.print_matrix()
        check = [
            ['1.000', '1.000', '1.000', '1.000', '1.000', '1.000'],
            ['1.000', '0.500', '0.000', '0.000', '-0.500', '-1.000'],
            ['1.000', '0.000', '1.000', '-1.000', '0.000', '1.000'],
            ['1.000', '0.000', '-1.000', '0.000', '0.000', '1.000'],
            ['1.000', '-0.500', '0.000', '0.000', '0.500', '-1.000'],
            ['1.000', '-1.000', '1.000', '1.000', '-1.000', '1.000']]
        assert(xis == check)

    def _test_octahedron(self):
        Octahedron.print_matrix()

    def _test_fct_cell(self):
        fct = FaceCenterTetragonal(2)
        gxi_oct = fct.get_octa_xi_glob()
        gxi_tet = fct.get_tetra_xi_glob()
        print(gxi_tet)

    def test_montecalro(self):
        cell = FaceCenterTetragonal(10, conc=0.95)
        ecis = np.array(
            [-11.370633, -611.080423, 688.864922, 271.55793, -10.2773419,
             -59.7960943, 0, 195.368747, 209.502264, 108.438473, 182.923418,
             -663.075602, -786.354929, -164.142113, 68.6509158, 648.47156, 0,
             -44.4990322]) * 1/1000.
        ecis[1] = -1300 * 1/1000.
        moncal = MonteCarlo(ecis, cell, 50)
        moncal._iterationTO(40)
        cell.make_poscar('POSCAR')

    def _test_energy_TO(self):
        _ecis = [688.864922, 108.438473, -663.075602, 182.923418,
                -44.4990322, -10.2773419, -11.3706330, -59.7960943,
                195.368747, 209.502264, -786.354929, 271.557930,
                -611.080423, -164.142113, 648.471560, 68.6509158]
        order = [2, 9, 11, 10, 17, 4, 0, 5, 7, 8, 12, 3, 1, 13, 15, 14]
        ecis = []
        for i in range(18):
            if i in order:
                pos = order.index(i)
                ecis.append(_ecis[pos])
            else:
                ecis.append(0)

        print('L10')
        fct = FaceCenterTetragonal(3, 'L10')
        print((fct.get_xi_glob() * ecis).sum())
        print()
        print('L12 A3B1')
        fct = FaceCenterTetragonal(3, 'L12_A')
        print((fct.get_xi_glob() * ecis).sum())
        print()
        print('L12 A1B3')
        fct = FaceCenterTetragonal(3, 'L12_B')
        print((fct.get_xi_glob() * ecis).sum())

        print(ecis)

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
