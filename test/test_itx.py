#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for itx.py"""
import os
import glob
import shutil
import unittest

import itx

TEST_PATH = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
             "99_python/01_testRun/")

def main():
    unittest.main()

class TestWave(unittest.TestCase):
    """
    itxを作成
    """
    path = os.path.join(TEST_PATH, 'igor')
    def test_to_str(self):
        wave = itx.Wave('wave0', [1, 2, 3])
        print(wave.to_itx())
        print(type(wave) == itx.Wave)

    def test_pref(self):
        mu = r"\\F'Symbol'm\\F'Times New Roman'"
        wavex = itx.Wave('wavex', [1, 2, 3], r"\\Z20c/a")
        wave0 = itx.Wave('wave0', [2, 3, 4],
                         r"\\Z20MAE ({0}eV/atom)".format(mu))
        wave2 = itx.Wave('wave2', [1, 3, 5], r"\\Z20Enthalpy (kJ/atom)")
        wave1 = itx.Wave('wave1', [3, 1, 0],
                         r"\\Z20Mag. ({0}\\M\\Z12B\\Z20/atom)".format(mu))

        waves = [wavex, wave0, wave1, wave2]
        pref = itx.Produce('test_plot', waves)
        pref.vertical3(waves)
        print(pref.to_itx())
        fname = os.path.join(self.path, 'test.itx')
        with open(fname, 'w') as wfile:
            wfile.write(pref.to_itx())

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
