#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""volume依存性"""
import os
import re
import unittest
import pylab

TEST_PATH = "/Users/enoki/Researches/Analysis/Codes/01_testRun/"

def main():
    unittest.main()

class iCVM_str(object):
    """
    iCVM の **.str file を取り扱う
    四元系
    """
    def __init__(self):
        pass

    @classmethod
    def from_file(cls, fname):
        """
        file から object を作成
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        print(len(lines))
        meta_head = re.compile(r"(.*)DISORD=(.*)")
        j = 0
        for i in range(len(lines)):
            if meta_head.match(lines[i]):
                formula = meta_head.match(lines[i]).group(1).split()[0].split("*")[1]
                meta = re.compile(r"([A-Z][0-9]+)"*4)
                comp = [int(meta.match(formula).group(x)[1:])
                        for x in range(1,5)]
                cell = lines[i+1:i+4+sum(comp)]
                fraction = comp[0]/(comp[0] + comp[2] + comp[3])
                if fraction < 0.26 and fraction > 0:
                    j += 1
        print(j)

class testParse(unittest.TestCase):
    path = os.path.join(TEST_PATH, 'icvm', 'str')

    def test_iCVM_str(self):
        """
        iCVM_str class の test
        """
        fname = os.path.join(self.path, "bcci.str")
        iCVM_str.from_file(fname)

if __name__ == '__main__':
    main()
