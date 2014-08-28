#!/usr/bin/python
# -*- coding: utf-8 -*-
"""test"""
import os
import glob
import pylab
import shutil
import unittest
import uspex
import convex_hull
from commopy import Bash


TEST_PATH = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
             "99_python/01_testRun/")

def main():
    unittest.main()

class TestUspexPOSCARS(unittest.TestCase):
    """
    POSCARSを展開するテスト
    """
    path = os.path.join(TEST_PATH, 'uspex', 'POSCAR')
    def test_split_structure(self):
        """
        構造毎に別けて出力
        Fe-Bの二元系計算が例
        Feが0.94<x<1の範囲の組成のみを抽出
        """
        clean_prev_dir(self.path, 'EA*')

        fname = os.path.join(self.path, 'extended_convex_hull_POSCARS')
        poscars = uspex.POSCARS(fname)

        poscars.restrict_fractions(0, 0.94)
        poscars.restrict_fractions(1, 0.0)
        for poscar in poscars.poscars:
            print(poscar['object'].poscar_title)
            print(poscar['object'].get_atom_fractions())
        print(len(poscars.poscars))
        print(poscars.poscars[0]['object'].poscar_title)
        print(poscars.poscars[-1]['object'].poscar_title)
        poscars.remove_elements(1)
        #poscars.expand_poscars(['Fe', 'B'])

        poscars.expand_poscars(['Fe'])


class TestOrigin(unittest.TestCase):
    path = os.path.join(TEST_PATH, 'uspex', 'origin')

    def test_draw_family_tree(self):
        path_orig = os.path.join(self.path, 'origin')
        origin = uspex.Origin(path_orig, 1866)
        origin.draw_family_tree(self.path)

class TestAuxiliary(unittest.TestCase):
    path2 = os.path.join(TEST_PATH, 'uspex', 'binary')
    path3 = os.path.join(TEST_PATH, 'uspex', 'ternary')

    def test_binary(self):
        aux = uspex.Auxiliary.from_file(self.path2)
        aux.binary()

    def test_ternary(self):
        out = uspex.Output(os.path.join(self.path3, 'OUTPUT.txt'))
        aux = uspex.Auxiliary.from_file(self.path3)
        aux.ternary()
        initial_base, not_bases, meta_stables = aux.separate_bases()
        ax = convex_hull.draw_convex_hull(initial_base, not_bases,
                                          out.elements, [-60, 5])
        pylab.show()

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
