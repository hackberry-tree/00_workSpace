#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
cvm imput file を修正したりする

e.g.
    C1, D1のラベルを入れ替える:
        cvm_imput.py bcci.str --exchange C1 D1

    元素を追加:
        cvm_imput.py bcci.str --add_site input
        >> input の中身
            primitive_z = 3
            origin = 0
            trans = -1/2 0 1/2
            label = B0
            insert_idx = 2

    A の濃度が 0.25 より小さいもののみを収集
        cvm_imput.py bcci.str --reject 0.25

"""
from __future__ import division

import argparse
from fractions import Fraction
from parse_cvm import CVMStrs, CVMEnergies

def main():
    """
    use argparse
    """
    parser = argparse.ArgumentParser(description='cvm parser')
    parser.add_argument('path', type=str, nargs=1)

    # exchange C and D
    parser.add_argument('--exchange', dest='exchange', type=str, nargs=2,
                        action='store', default=None)

    parser.add_argument('--add_site', dest='add_site', type=str, nargs=1,
                        action='store', default=None)

    parser.add_argument('--reject', dest='reject', type=float, nargs=1,
                        action='store', default=None)

    args = parser.parse_args()
    opt = {'path': args.path[0]}
    if args.exchange:
        exe = exchange_elements
        opt.update({'elem1': args.exchange[0], 'elem2': args.exchange[1]})

    if args.add_site:
        exe = add_site
        opt.update({'input_file': args.add_site[0]})

    if args.reject:
        exe = reject
        opt.update({'limit': args.reject[0]})

    exe(**opt) #pylint: disable=W0142

def exchange_elements(path, elem1, elem2):
    """
    原子ラベルを交換したファイルを作成
    """
    strs = CVMStrs.from_str_file(path)
    enes = CVMEnergies.from_energies_file("energies.txt")
    strs.exchange_elements(elem1, elem2)
    enes.exchange_elements(elem1, elem2)
    fname = path + "_rev"
    strs.make_file(fname)
    enes.make_file("energies.txt_rev")

def add_site(path, input_file):
    """
    サイトを追加
    """
    with open(input_file, 'r') as rfile:
        lines = rfile.readlines()
    args = {}
    for line in lines:
        key = line.split("=")[0].split()[0]
        tmp = line.split("=")[1].split()
        if len(tmp) == 1:
            try:
                value = int(tmp[0])
            except ValueError:
                value = tmp[0]
        else:
            value = [Fraction(x) for x in tmp]
        args.update({key: value})
    strs = CVMStrs.from_str_file(path)
    strs.add_site(**args) #pylint: disable=W0142
    fname = path + "_added"
    strs.make_file(fname)

def reject(path, limit):
    """
    条件で構造を絞る
    """
    strs = CVMStrs.from_str_file(path)
    tmpstr = []
    for structure in strs:
        a = structure.label.num_atoms[0] / (structure.label.num_atoms[0] +
                                            structure.label.num_atoms[2] +
                                            structure.label.num_atoms[3])
        if a != 0 and a < limit:
            tmpstr.append(structure)
    outstr = CVMStrs(tmpstr)
    fname = path + "_rejcted"
    outstr.make_file(fname)

        # print(a)


if __name__ == '__main__':
    main()
