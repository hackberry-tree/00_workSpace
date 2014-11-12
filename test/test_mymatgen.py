#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for mymatgen.py"""
import os
import glob
import pickle
import shutil
import unittest

import mymatgen
from makeSeries import Combinatorial

TEST_PATH = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
             "99_python/01_testRun/")

def main():
    unittest.main()

class TestMyMat(unittest.TestCase):
    """
    test_endmemb: end_membの組成が正しいかcheck
    """
    #path = os.path.join(TEST_PATH, 'igor')
    source_path = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
                   "99_python/00_workSpace/sorce/pymatgen/ternary_convex")
    def test_endmemb(self):
        elem_list = [['Fe'],
                     ['Al', 'Co', 'Cr', 'Cu', 'Ga',
                      'Ge', 'Mn', 'Nb', 'Ni', 'Pd', 'Pt',
                      'Rh', 'Ru', 'Sb', 'Si', 'Sn', 'Ti',
                      'V', 'Zn', 'Zr'],
                     ['C', 'B', 'N']]
        elem_list = [['Fe'], ['Ni', 'Ge', 'Ga', 'Co'], ['Si', 'Al']]
        elem_list = [['Fe'], ['Nb'], ['B', 'C']]
        combi = Combinatorial(*elem_list)
        comp_list = [x['elements'] for x in combi.compositions]
        print(len(comp_list))
        for comp in comp_list:
            pickle_name = "".join(comp) + ".pickle"
            print(pickle_name)
            pickle_path = os.path.join(self.source_path, pickle_name)
            if glob.glob(pickle_path):
                with open(pickle_path, 'rb') as rbfile:
                    end_memb = pickle.load(rbfile)[0]
            else:
                mp_data = mymatgen.Ternary(comp)
                end_memb = mp_data.to_convex_hull()[0]
            end_comp = [[0, 0], [1.0, 0], [0, 1.0]]
            self.assertTrue(end_memb[0][0:2] in end_comp)
            self.assertTrue(end_memb[1][0:2] in end_comp)
            self.assertTrue(end_memb[2][0:2] in end_comp)

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
