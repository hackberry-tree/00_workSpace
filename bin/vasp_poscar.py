#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
POSCARを取り扱う
"""
from __future__ import division

import argparse
from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifWriter
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

    # make cif file
    parser.add_argument('--cif', dest='run', const=cif,
                        action='store_const', default=None)

    # convert to refined POSCAR
    parser.add_argument('--refined', dest='run', const=refined,
                        action='store_const', default=None)

    # print spg
    parser.add_argument('--print_spg', dest='run', const=print_spg,
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
    std = finder.get_conventional_standard_structure()
    dstpos = Poscar(std)
    dst = 'POSCAR_std'
    Cabinet.reserve_file(dst)
    dstpos.write_file(dst)

def cif(src='POSCAR'):
    """
    cifファイルを作成
    """
    srcpos = Poscar.from_file(src)
    finder = SpacegroupAnalyzer(srcpos.structure)
    std = finder.get_conventional_standard_structure()
    cif = CifWriter(std, find_spacegroup=True, symprec=0.1)
    cif.write_file('poscar.cif')

def refined(src='POSCAR'):
    srcpos = Poscar.from_file(src)
    finder = SpacegroupAnalyzer(srcpos.structure,
                                symprec=5e-2, angle_tolerance=8)
    std = finder.get_refined_structure()
    dstpos = Poscar(std)
    dst = 'POSCAR_refined'
    Cabinet.reserve_file(dst)
    dstpos.write_file(dst)

def print_spg(src='POSCAR'):
    srcpos = Poscar.from_file(src)
    finder = SpacegroupAnalyzer(srcpos.structure,
                                symprec=5e-2, angle_tolerance=8)
    spg = finder.get_spacegroup_symbol()
    spg_num = finder.get_spacegroup_number()
    print(spg)
    print(spg_num)

if __name__ == '__main__':
    main()
