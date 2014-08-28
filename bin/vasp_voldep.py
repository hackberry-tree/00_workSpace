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
    parser.add_argument('--5to10', dest='run', )
    parser.add_argument('--run', dest='run', const=from5to10)
    args = parser.parse_args()
    args.run()


def main():
    from5to10()

def frange():
    volume_list = Array.frange_stp(90, 140, 5)
    series = series_vasp.Produce('POSCAR', 'vol')
    series.set_volume(volume_list)
    series.make_files()
    series.append_list_run('run_volume.sh')

def from5to10():
    poscar = vaspy.Poscar()
    volume = poscar.get_cell_volume()
    volume_list = [volume*(0.95+i*0.15/10) for i in range(0, 11)]
    series = series_vasp.Produce('POSCAR', 'voldep')
    series.set_volume(volume_list)
    series.make_files()
    series.append_list_run('run_volume.sh')

if __name__ == '__main__':
    main()
