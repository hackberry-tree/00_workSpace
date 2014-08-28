#!/usr/bin/python
# -*- coding: utf-8 -*-
"""test"""
import os
import shutil
import unittest
from commopy import Bash, Table
import numpy as np



TEST_PATH = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
             "99_python/01_testRun/")

def main():
    unittest.main()


class testTable(unittest.TestCase):
    path = os.path.join(TEST_PATH)

    def test_table(self):
        data = Table([{},{}])
        data.update('test', [{'Fe':1, 'Co':2},{'Fe':12, 'Co':30}])
        print(data)

def clean_prev(path, files):
    """
    既存filesを消去
    """
    trush_list = Bash.find_files(path, files)
    for trush in trush_list:
        fname = os.path.join(path, trush)
        os.remove(fname)
        print("{0} is removed.".format(fname))


def clean_prev_dir(path, dirc):
    """
    既存dirctoryを消去
    """
    trush_list = Bash.find_files(path, dirc)
    for trush in trush_list:
        fname = os.path.join(path, trush)
        shutil.rmtree(fname)
        print("{0} is removed.".format(fname))


if __name__ == '__main__':
    main()
