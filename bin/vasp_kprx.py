#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
kpoints依存性等間隔にふれないのであまり良くない
"""
import series_vasp
from commopy import Bash, Array


def main():
    Bash.mkdir('kp')
    kp_list = Array.frange_stp(0.08, 0.24, 0.02)
    #kp_list = Array.frange_stp(0.04, 0.07, 0.02)
    series = series_vasp.Produce('POSCAR', 'kp')
    series.set_kp_relax(kp_list)
    series.make_files()
    series.make_list_run('run.sh')


if __name__ == '__main__':
    main()
