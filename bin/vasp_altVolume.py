#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
組成に応じて体積の初期値を修正する
i-s の構造ファイル用

① 組成を POSCAR から読む
② 体積計算
③ 変更
"""

from __future__ import division

import os
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
    use argparse
    """
    parser = argparse.ArgumentParser(description='convert POSCAR')
    parser.add_argument('--check', dest='check', type=str, nargs=1,
                        action='store', default=None)

    parser.add_argument('--alt', dest='run', const=alt_volume,
                        action='store_const', default=None)

    args = parser.parse_args()
    if args.check:
        # print(args.check[0])
        check(args.check[0])
        return

    args.run()

def get_elements(src):
    srcpos = Poscar.from_file(src)
    elements = [x.symbol for x in srcpos.structure.species]
    elem_counter = Counter(elements)
    return elem_counter

def calc_volume(elem_counter):
    Fe = elem_counter['Fe']
    Cr = elem_counter['Cr']
    subs = Fe + Cr

    C = elem_counter['C']

    Fe_C = 22.863 * (1 - math.exp(-0.52018 * C / (subs)))
    Cr_C = 26.839 * (1 - math.exp(-0.49806 * C / (subs)))

    vol = []
    vol.append((12.518 + Fe_C) * Fe)
    vol.append((12.870 + Cr_C) * Cr)
    return sum(vol)

def alt_volume(src='POSCAR'):
    srcpos = Poscar.from_file(src)
    elem = get_elements(src)
    ideal_v = calc_volume(elem)
    srcpos.structure.scale_lattice(ideal_v)
    dst = 'POSCAR_newvol'
    srcpos.write_file(dst)



def check(path):
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    lines = ""
    for dirc in dir_list:
        elem = get_elements(os.path.join(dirc, 'POSCAR'))
        lines += str(elem['Fe']) + "\t" + str(elem['Cr']) + "\t" + str(elem['C'])
        lines += "\t" + str(calc_volume(elem)) + "\n"
    print(lines)





if __name__ == '__main__':
    main()
