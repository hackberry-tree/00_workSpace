#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""cova依存性"""
import series_vasp
from commopy import Bash, Array


def main():
    Bash.mkdir('cova')
    cova_list = Array.frange_stp(1.10, 1.35, 0.025)
    cova_list.append(1.00)
    cova_list.append(1.05)
    series = series_vasp.Produce('POSCAR', 'cova')
    series.set_cova(cova_list)
    series.make_files()
    series.append_list_run('run_cova.sh')

if __name__ == '__main__':
    main()
