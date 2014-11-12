#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""cova依存性"""
import series_vasp
from commopy import Bash, Array


def main():
    Bash.mkdir('cova')
    #cova_list = Array.frange_stp(0.95/(2**0.5), 1.45/(2**0.5), 0.05/(2**0.5))

    # from bcc to fcc
    #cova_list = Array.frange_stp(0.90, 1.60, 0.05)

    # from bcc to fcc for D03
    rate = 2 ** 0.5
    rate = 1. / (2 ** 0.5)
    rate = 1.

    cova_list = Array.frange_stp(0.90*rate, 1.60*rate, 0.05*rate)
    #cova_list = [1.65, 1.7, 1.75, 1.8]

    #cova_list.append(1.00)
    #cova_list.append(1.05)
    series = series_vasp.Produce('POSCAR', 'cova')
    series.set_cova(cova_list)
    series.make_files()
    series.append_list_run('run_cova.sh')

if __name__ == '__main__':
    main()
