#!/usr/bin/python
# -*- coding: utf-8 -*-
"""test"""
import os
import unittest
import judgeRX


def main():
    unittest.main()


class VaspyIncar(unittest.TestCase):
    path = ('/Users/enoki/Documents/01_ResearchData/Calculations'
            '/99_python/01_testRun/judge_relax')
    os.chdir(path)
    def test_judgeOsz(self):
        print("test_judgeOsz")
        judge = judgeRX.judge_oszicar('1finish')
        self.assertTrue(judge)
        judge = judgeRX.judge_oszicar('morethanNSW')
        self.assertFalse(judge)
        judge = judgeRX.judge_oszicar('smalldiff')
        self.assertTrue(judge)


if __name__ == '__main__':
    main()
