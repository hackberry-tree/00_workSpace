#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for parse_atat"""


import os
import glob
import shutil
import unittest

import parse_atat

__date__ = "Sep 3 2014"

TEST_PATH = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
             "99_python/01_testRun/")

def main():
    unittest.main()

class Test(unittest.TestCase):
    """
    テスト
    """
    PATH = os.path.join(TEST_PATH, 'atat')

    def _test_strout(self):
        """
        StrOutのテスト
        """
        src = os.path.join(self.PATH, 'Fe-Ni', '10487')

        #strout = parse_atat.StrOut.from_file(src)
        #print(strout.structure)
        #print(dir(strout.prim_cif))
        #strout.prim_cif
        #dst = os.path.join(self.PATH, 'Fe-Ni', '10487', 'str.cif')
        #strout.prim_cif(dst)

    def test_analysis(self):
        dirc = os.path.join(self.PATH, 'Fe-Ni')
        analysis = parse_atat.Analysis.from_dirc(dirc)
        print(analysis.to_tex_form())




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
