#!/usr/bin/python
# -*- coding: utf-8 -*-
"""test"""
import os
import pylab
import shutil
import unittest
import convex_hull
from commopy import Bash


TEST_PATH = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
             "99_python/01_testRun/")

def main():
    unittest.main()


class TestConvex(unittest.TestCase):
    path = os.path.join(TEST_PATH, 'convex_hull', 'FeNiSi')

    def test_draw_convex_hull(self):
        data = pylab.loadtxt(self.path, comments='#')
        initial_base = [list(x) for x in data[0:3, [0, 2, 3]]]
        not_base = [list(x) for x in data[3:, [0, 2, 3]]]
        convex_hull.draw_convex_hull(initial_base, not_base, ['Fe', 'Ni', 'Si'], [-60, 5])
        print(not_base)
        pylab.show()

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
