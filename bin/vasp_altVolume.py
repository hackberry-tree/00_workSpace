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

def main():
    """
    use argparse
    """
    parser = argparse.ArgumentParser(description='convert POSCAR')
    parser.add_argument('--check', dest='check', type=str, nargs=2,
                        action='store', default=None)

    parser.add_argument('--alt', dest='calc', type=str, nargs=1,
                        action='store', default=None)

    args = parser.parse_args()
    if args.check:
        # print(args.check[0])
        if args.check[1] == 'bcci':
            check(args.check[0], bcci)
        if args.check[1] == 'fcc':
            check(args.check[0], fcc)
        return

    if args.calc:
        if args.calc[0] == 'bcci':
            alt_volume(bcci)
        if args.calc[0] == 'fcc':
            alt_volume(fcc)

def get_elements(src):
    """
    組成比を得るための method
    """
    srcpos = Poscar.from_file(src)
    elements = [x.symbol for x in srcpos.structure.species]
    elem_counter = Counter(elements)
    return elem_counter

def bcci(elem_counter):
    """
    体積を計算する
    bcci 用 (Fe, Cr, Ti, C)
    """
    Fe = elem_counter['Fe']
    Cr = elem_counter['Cr']
    Ti = elem_counter['Ti']
    subs = Fe + Cr + Ti

    C = elem_counter['C']

    Fe_C = 22.863 * (1 - math.exp(-0.52018 * C / (subs)))
    Cr_C = 26.839 * (1 - math.exp(-0.49806 * C / (subs)))
    Ti_C = 35.809 * (1 - math.exp(-0.30445 * C / (subs)))

    vol = 0
    vol += (12.518 + Fe_C) * Fe
    vol += (12.870 + Cr_C) * Cr
    vol += (16.949 + Ti_C) * Ti
    return vol

def fcc(elem_counter):
    """
    体積を計算する
    fcc 用 (Al, Cu)
    """
    Al = elem_counter['Al']
    Cu = elem_counter['Cu']
    vol = 0
    vol += 16.154 * Al
    vol += 11.4594 * Cu
    return vol

def alt_volume(func, src='POSCAR'):
    """
    体積を変更した POSCAR を作成
    """
    srcpos = Poscar.from_file(src)
    elem = get_elements(src)
    ideal_v = func(elem)
    srcpos.structure.scale_lattice(ideal_v)
    dst = 'POSCAR_newvol'
    srcpos.write_file(dst)

def check(path, func):
    """
    動作チェック用
    """
    print(path)
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    print(dir_list)
    lines = ""
    for dirc in dir_list:
        elem = get_elements(os.path.join(dirc, 'POSCAR'))
        # lines += str(elem['Fe']) + "\t" + str(elem['Cr']) + "\t" + str(elem['C'])
        # lines += str(elem['Fe']) + "\t" + str(elem['Ti']) + "\t" + str(elem['C'])
        lines += str(elem['Al']) + "\t" + str(elem['Cu'])
        lines += "\t" + str(func(elem)) + "\n"
    print(lines)





if __name__ == '__main__':
    main()
