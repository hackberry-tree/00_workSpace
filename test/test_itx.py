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
    POSCARSを展開するテスト
    """
    def test_to_str(self):
        wave = itx.Wave('wave1', [1, 2, 3])
        print(wave)

    def test_pref(self):
        wave = itx.Wave('wave1', [1, 2, 3])
        waves = [wave, wave, wave ,wave]
        pref = itx.Preferences()
        pref.vertical3(waves)
        print(pref.lines)


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
