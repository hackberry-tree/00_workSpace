#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""volume依存性"""
import glob
import collect_vasp
from commopy import Bash, Array

def main():
    path_list = glob.glob('volume_*')
    data = collect_vasp.Energy(path_list, 'POSCAR', 'OSZICAR')
    data.data.sort(key=lambda x: x['volume'])
    data.output_keys = ['volume', 'energy', 'mag']
    print(data)

if __name__ == '__main__':
    main()
