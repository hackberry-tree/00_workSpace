#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
POSCARSを展開
"""
import os
import uspex
from commopy import Cabinet

def main():
    """main"""
    #for_FeB()
    for_FeB_witoutB()

def for_FeB(): #pylint: disable=C0103
    """
    FeBのUSPEXの計算結果(POSCARS)を展開する
    Feが50-100percentの範囲
    """
    fname = os.path.join('.', 'extended_convex_hull_POSCARS')
    poscars = uspex.POSCARS(fname)

    poscars.restrict_fractions(0, 0.5)
    poscars.restrict_fractions(1, 0.0)
    poscars.expand_poscars(['Fe', 'B'])
    path_list = [x['ID'] for x in poscars.poscars]
    Cabinet.make_list_run(path_list, 'run.sh')

def for_FeB_witoutB(): #pylint: disable=C0103
    """
    for_FeBのBのサイトを0にして計算
    """
    fname = os.path.join('.', 'extended_convex_hull_POSCARS')
    poscars = uspex.POSCARS(fname)

    poscars.restrict_fractions(0, 0.5)
    poscars.restrict_fractions(1, 0.0)
    poscars.remove_elements(1)
    poscars.expand_poscars(['Fe'])
    path_list = [x['ID'] for x in poscars.poscars]
    Cabinet.make_list_run(path_list, 'run.sh')


if __name__ == '__main__':
    main()
