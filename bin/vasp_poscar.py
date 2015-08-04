#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
POSCARを取り扱う
"""
from __future__ import division

import os
import glob
import argparse
from collections import Counter
from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_analyzer import RelaxationAnalyzer

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

    parser.add_argument('--fname', dest='fname', type=str, nargs=1,
                        action='store', default=['POSCAR'])
    # five
    parser.add_argument('--five', dest='run', const=five,
                        action='store_const', default=None)
    # checkrelax
    parser.add_argument('--chrx', dest='run', const=checkrelax,
                        action='store_const', default=None)
    # volume
    parser.add_argument('--volume', dest='run', const=volume,
                        action='store_const', default=None)

    parser.add_argument('--newv', dest='newv', type=float, nargs=1,
                        action='store', default=[None])

    args = parser.parse_args()
    print(args.fname[0])

    if args.run == volume:
        args.run(args.fname[0], args.newv[0])
        return

    args.run(args.fname[0])

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

def five(src='POSCAR'):
    """
    ver5 に変換
    """
    srcpos = Poscar.from_file(src)
    dst = 'POSCAR_five'
    srcpos.write_file(dst)

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
    """
    refined poscar を 作成する
    """
    srcpos = Poscar.from_file(src)
    finder = SpacegroupAnalyzer(srcpos.structure,
                                symprec=5e-1, angle_tolerance=8)
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

def checkrelax_single(path, src_ini='posfinal', src_fin='posfinal3'):
    dirc = path
    initial = Poscar.from_file(
        os.path.join(dirc, src_ini)).structure.as_dict()['lattice']
    final = Poscar.from_file(
        os.path.join(dirc, src_fin)).structure.as_dict()['lattice']
    length_a = final[u'a']/initial[u'a']
    length_b = final[u'b']/initial[u'b']
    length_c = final[u'c']/initial[u'c']
    delta_length = ((length_b / length_a - 1) ** 2 +
                    (length_c / length_a - 1) ** 2) ** 0.5
    # print(delta_length)
    angle_a = final['alpha']/initial['alpha'] - 1
    angle_b = final['beta']/initial['beta'] - 1
    angle_c = final['gamma']/initial['gamma'] - 1
    delta_angle = (angle_a ** 2 + angle_b ** 2 + angle_c ** 2) ** 0.5
    # print(delta_angle)
    return delta_length, delta_angle

def checkrelax(path):
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    lines = ""
    for dirc in dir_list:
        dist = checkrelax_single(dirc)
        lines += os.path.basename(dirc) + "\t"
        lines += str(dist[0]) + "\t" + str(dist[1]) + "\n"
    print(lines)

def volume(src, new_volume=None):
    """
    volume を変更
    単位は A^3 per atom
    """
    srcpos = Poscar.from_file(src)
    pre_volume = srcpos.structure.volume
    print('present volume')
    print(pre_volume)
    if new_volume:
        srcpos.structure.scale_lattice(new_volume)
        print('')
        print('new volume (POSCAR_newvol)')
        print(new_volume)
    dst = 'POSCAR_newvol'
    srcpos.write_file(dst)

def get_elements(src):
    srcpos = Poscar.from_file(src)
    elements = [x.symbol for x in srcpos.structure.species]
    elem_counter = Counter(elements)
    return elem_counter




if __name__ == '__main__':
    main()
