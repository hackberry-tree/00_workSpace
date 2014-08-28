#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fix no member in numpy
"""
import os
import re
import math
import copy
import numpy as np
import solid
from commopy import Cabinet, Vector, Bash



def main():
    modname = 'numpy.random'
    fname = "/opt/local/Library/Frameworks/Python.framework/Versions/3.3/lib/python3.3/site-packages/numpy/random/__init__.py"
    fix(modname, fname)

def fix(modname, fname):
    """
    modnameはパッケジされている場合に必要
    """
    lines = Cabinet.read_file(fname)
    key = r"\s*from\s+(.*)\s+import\s+\*"
    meta = re.compile(key)
    star_list = [meta.match(x).group(1) for x in lines if meta.match(x)]
    for module in star_list:
        exec('import {0}{1}'.format(modname, module))
        try:
            a = eval('{0}{1}.__all__'.format(modname, module))
        except AttributeError:
            b = eval('dir({0}{1})'.format(modname, module))
            a = [x for x in b if x[0] != '_']
        print('from {0} import {1}'.format(module, ', '.join(a)))
        print()



if __name__ == '__main__':
    main()
