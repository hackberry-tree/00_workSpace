#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MAGMOM を ferri 磁性 に修正
"""

from __future__ import division

import os
import re
import glob
import math
import argparse
from collections import Counter
from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_analyzer import RelaxationAnalyzer

from commopy import Cabinet

def main():
    """
    """
    get_elements()

def get_elements():
    with open("POSCAR", 'r') as rfile:
        lines = rfile.readlines()
    elements = lines[5].split()
    num_elem = lines[6].split()
    spin = {"Fe": 2.5, "Cr": -1, "C": 0, "Ti": 0, "N": 0}
    mag = [i + "*" + str(spin[e]) for e, i in zip(elements, num_elem)]
    magmom = "MAGMOM = " + "  ".join(mag) + "\n"

    with open("INCAR", 'r') as rfile:
        lines = rfile.readlines()
    meta = re.compile(r".*MAGMOM.*")
    i = 0
    while not meta.match(lines[i]):
        i += 1
    lines[i] = magmom
    with open("INCAR_ferri", 'w') as wfile:
        wfile.write("".join(lines))

if __name__ == '__main__':
    main()
