#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""組成依存性"""
import series_vasp
from commopy import Bash
from makeSeries import Combinatorial

def main():
    """main"""
    for_perovskite()

def for_ferrite():
    """
    ferrite系のcombinatorial
    """
    Bash.mkdir('combi')
    elem_list = ['Ni', 'Mn', 'Cu', 'Zn', 'Cr']
    comp_list = [['O', 'Fe', x] for x in elem_list]
    series = series_vasp.Produce('POSCAR', 'combi')
    series.set_elements(comp_list)
    series.make_files()
    series.append_list_run('run_combi.sh')


def for_perovskite():
    """
    perovskite系のcombinatorial
    """
    Bash.mkdir('combi')
    elem_list = [['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Ga',
                  'Ge', 'Mn', 'Nb', 'Ni', 'Pd', 'Pt',
                  'Rh', 'Ru', 'Sb', 'Si', 'Sn', 'Ti',
                  'V', 'Zn', 'Zr'], ['C', 'B', 'N']]
    combi = Combinatorial(*elem_list)
    comp_list = [['Fe', x['elements'][0], x['elements'][1]]
                 for x in combi.compositions]
    series = series_vasp.Produce('POSCAR', 'combi')
    series.set_elements(comp_list)
    series.make_files()
    series.append_list_run('run_combi.sh')


if __name__ == '__main__':
    main()
