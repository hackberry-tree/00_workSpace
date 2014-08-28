#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import math
import re
import copy
import itertools
import subprocess
import os
from trial_run import TrialRun

def main():
    a = test()
    dirc = '.'
    process = subprocess.Popen(['ls', '-a'], stdout=subprocess.PIPE)
    out, err = process.communicate()
    print(out)
    #Trial = TrialRun(testRun)
    #Trial.whetherExist(os.path.join(dirc, 'POSCAR'))

def testRun():
    print(test)
    Prep = test
    Prep.pprint
    #print Prep.defaultInput

class test:
    def __init__(self):
        print('OK')

    def pprint(self):
        print('err')

class PrepInputs:
    def __init__(self):
        self.sorce_dir = os.path.dirname(os.path.abspath(__file__))
        print(self.sorce_dir)
        self.defaultInput = readFile(os.path.join(self.sorce_dir, 'originalAKAI', 'input'))

    def hoge(self):
        pass

def readFile(filename):
    f = open(filename,'r')
    lines = f.readlines()
    f.close
    return lines

def writeFile(fname, lout):
    out = open(fname,'w')
    for line in lout:
        out.write(line)
    out.close

if __name__ == '__main__':
    main()

