#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for parse_atat"""


import os
import glob
import timeit
import shutil
import unittest
import correlationfunction as cf
import parse_atat
# import pyximport; pyximport.install(pyimport=True)
# from mcc import *
from mc import *
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

    @staticmethod
    def _test_octahedron_xi():
        """
        TO の相関関数を print
        両者が一致するか確認
        """
        Octahedron.print_matrix()
        cf.Octahedron.print_matrix()

    @staticmethod
    def _test_tetrahedron_xi():
        """
        TT の相関関数を print
        両者が一致するか確認
        """
        Tetrahedron.print_matrix()
        cf.Tetrahedron.print_matrix()

    @staticmethod
    def _test_xi_TO():
        """
        tetrahedron octahedronの相関関数をチェック
        特に四面体からのnull, point, pair_nn, tri_nnと
        八面体のとが同じ値を取っているかtestする
        """
        cell = FaceCenterCubic(30, arrange='L12_B')
        assert all((cell.get_octa_xi_glob()[[0, 1, 2, 4]] -
                    cell.get_tetra_xi_glob()[: 4]) ** 2 < 1e-8)

    @staticmethod
    def _test_dxi():
        """
        スピンフリップに伴う相関関数の変化量を check する
        """
        cell = FaceCenterCubic(2, arrange='L12_A')
        print('concentration')
        print(cell.get_conc(1))
        print('gloval xi')
        print(cell.get_tetra_xi_glob())
        cuau = Tetrahedron(0.5)
        print("Delta Xi")
        xi = cell.get_tetra_dxi(0)
        print(xi.shape)
        print(cuau.ecis.shape)
        de = xi * cuau.ecis
        print(np.exp(de.sum(axis=1)))
        print(np.exp(1))

    @staticmethod
    def _test_orderparam():
        """
        order parameter 計算の test
        """
        cell = FaceCenterCubic(2, arrange='L12_A')
        cell = FaceCenterCubic(24, arrange='random')
        cuau = Tetrahedron(0.5).ecis
        moncal = MonteCarlo(cuau, cell, 100)
        #de = moncal.delta_e()
        #print(de)
        #print(moncal.transition_probability(de))
        print(moncal._iteration(10))
        print(cell.get_orderparam(1))
        print(cell.get_orderparam(2))
        print(cell.get_orderparam(3))
        print(cell.get_orderparam(4))
        print(cell.get_orderparam(6))
        print(cell.get_orderparam(11))

    def test_gs_TO(self):
        """
        TO近似の規則相のエネルギー比較 check
        """
        path = os.path.join(TEST_PATH, "montecalro/", "AlCu", "wien", "TO")
        cell = FaceCenterCubic(10, arrange='L10', conc=0.95)
        fname = os.path.join(path, 'POSCAR.pickle')
        #cell = FaceCenterCubic.load_cell(fname)

        for i in range(1):
            #dirc = os.path.join(path, 'data_{0}'.format(i))
            #os.makedirs(dirc)
            t = 500/(i+1)
            t = 300
            moncal = MonteCarlo(ecis, cell, t)
            #conc_ene = moncal._iterationTO_reserved_atoms(100)
            conc_ene = moncal._iterationTO(10)
            fname = os.path.join(path, 'POSCAR')
            cell.make_poscar(fname)
            stack_data.append([i, conc_ene[-1][0], conc_ene[-1][1]])
        cell.save_cell(fname)
        return

        lines = "\n".join([" ".join([str(x) for x in data])
                           for data in stack_data])
        fname = os.path.join(path, 'comp_ene.txt')
        with open(fname, 'w') as wfile:
            wfile.write(lines)
        cell.save_cell(fname)
        fname = os.path.join(path, 'POSCAR')
        cell.make_poscar(fname)

        return
        for i in range(1):
            dirc = os.path.join(path, 'data_{0}'.format(i))
            os.makedirs(dirc)

            cell = FaceCenterCubic(10, arrange='random')
            moncal = MonteCarlo(ecis, cell, 10)
            conc_ene = moncal._iterationTO(10)
            fname = os.path.join(dirc, 'POSCAR')
            cell.make_poscar(fname)
            stack_data.append([i, conc_ene[-1][0], conc_ene[-1][1]])
        lines = "\n".join([" ".join([str(x) for x in data])
                           for data in stack_data])
        fname = os.path.join(path, 'comp_ene.txt')
        with open(fname, 'w') as wfile:
            wfile.write(lines)


        # print(cell.get_orderparam(1))

    def _test_get_conc(self):
        cell = FaceCenterCubic(4, arrange='random')
        print(cell.get_conc(1))

    def _test_fct_print_coordinates(self):
        FaceCenterTetragonal.print_coordinates()

    def _test_fcc_print_coordinates(self):
        FaceCenterCubic.print_coordinates()

    def _test_fct(self):
        fct = FaceCenterTetragonal(4, arrange='random')

    def _test_fcc(self):
        fcc = FaceCenterCubic(4, arrange='random')
        print(fcc.OCTA == fcc._OCTA)


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
