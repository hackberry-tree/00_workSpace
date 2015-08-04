#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""volume依存性"""
import argparse
import series_vasp
import vaspy
from commopy import Bash, Array

def execute_parser():
    """
    parser
    """
    parser = argparse.ArgumentParser(description='make vol-dep calculation')
    parser.add_argument('--5to10', dest='run', const=from5to10,
                        action='store_const', default=None)

    parser.add_argument('--run', dest='runfile', type=str, nargs=1,
                        action='store', default=['run_volume.sh'])

    parser.add_argument('--recalc', dest='run', const=recalc_from_fit,
                        action='store_const', default=None)

    args = parser.parse_args()
    print(args.runfile[0])
    args.run(args.runfile[0])


def main():
    execute_parser()

def frange():
    volume_list = Array.frange_stp(90, 140, 5)
    series = series_vasp.Produce('POSCAR', 'vol')
    series.set_volume(volume_list)
    series.make_files()
    series.append_list_run('run_volume.sh')

def from5to10(runfile='run_volume.sh'):
    poscar = vaspy.Poscar()
    volume = poscar.get_cell_volume()
    volume_list = [volume*(0.95+i*0.15/10) for i in range(0, 11)]
    series = series_vasp.Produce('POSCAR', 'voldep')
    series.set_volume(volume_list)
    series.make_files()
    series.append_list_run(runfile)

def recalc_from_fit(runfile='run_volume.sh'):
    poscar = vaspy.Poscar()
    volume = poscar.get_cell_volume()

    fname='fit_results.dat'
    with open(fname, 'r') as rfile:
        lines = rfile.readlines()
    opt_volume = float(lines[3].split()[1])

    ratio = opt_volume * sum(poscar.num_atoms) / volume

    print(ratio)

    if ratio < 0.95:
        pt = (0.95 - ratio) // 0.015
        volume_list = [volume*(0.95-(i+1)*0.15/10) for i in range(int(pt+3))]
        series = series_vasp.Produce('POSCAR', 'voldep')
        series.set_volume(volume_list)
        series.make_files()
        series.append_list_run(runfile)

    if ratio > 1.10:
        pt = (ratio - 1.10) // 0.015
        volume_list = [volume*(1.10+(i+1)*0.15/10) for i in range(int(pt+3))]
        series = series_vasp.Produce('POSCAR', 'voldep')
        series.set_volume(volume_list)
        series.make_files()
        series.append_list_run(runfile)


if __name__ == '__main__':
    main()
