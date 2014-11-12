#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
POSCARを取り扱う
"""
from __future__ import division

import argparse
from pymatgen.io.vaspio import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from commopy import Cabinet

__date__ = "Oct 30 2014"


def main():
    """
    use argparse
    """
    parser = argparse.ArgumentParser(description='convert POSCAR')

    # convert to primitive
    parser.add_argument('--prim', dest='run', const=primitive,
                        action='store_const', default=None)

    # convert to standard
    parser.add_argument('--std', dest='run', const=standard,
                        action='store_const', default=None)

    args = parser.parse_args()
    args.run()

def primitive(src='POSCAR'):
    """
    primitiveに変換
    """
    srcpos = Poscar.from_file(src)
    finder = SpacegroupAnalyzer(srcpos.structure)
    prim = finder.get_primitive_standard_structure()
    dstpos = Poscar(prim)
    dst = 'POSCAR_prim'
    Cabinet.reserve_file(dst)
    dstpos.write_file(dst)

def standard(src='POSCAR'):
    """
    standardに変換
    """
    srcpos = Poscar.from_file(src)
    finder = SpacegroupAnalyzer(srcpos.structure)
    std = finder.get_primitive_standard_structure()
    dstpos = Poscar(std)
    dst = 'POSCAR_std'
    Cabinet.reserve_file(dst)
    dstpos.write_file(dst)


if __name__ == '__main__':
    main()
