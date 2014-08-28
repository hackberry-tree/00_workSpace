#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
mae計算
"""
import os
import sys
import glob
import collect_vasp

def main():
    """main"""
    if len(sys.argv) == 1:
        get_two_direction('.')
    elif sys.argv[1] == '--cova':
        get_cova()

def get_mae(base, alternative):
    """get mae"""
    path = ['.']
    energy = collect_vasp.Energy(path)
    energy.set_mae(base, alternative)
    energy.output_keys = ['mae']
    print(energy)

def get_two_direction(path='.'):
    """100-001 & 110-001"""
    print("100 - 001")
    soc001 = os.path.join(path, 'OSZICAR_soc001')
    soc100 = os.path.join(path, 'OSZICAR_soc100')
    soc110 = os.path.join(path, 'OSZICAR_soc110')
    get_mae(soc001, soc100)
    print("110 - 001")
    get_mae(soc001, soc110)

def get_cova():
    """
    cova dir中のデータを収集
    """
    dir_list = glob.glob("cova*/cova_*/OSZICAR_soc110")
    dir_list = [os.path.dirname(x) for x in dir_list]

    energy = collect_vasp.Energy(dir_list, output='OSZICAR_presoc_nc')
    energy.set_enthalpy()
    energy.set_mae('OSZICAR_soc001', 'OSZICAR_soc100')
    energy.data.sort(key=lambda x: x['c/a'])
    energy.output_keys = ['c/a', 'enthalpy', 'mae', 'mag']
    print(energy)


if __name__ == '__main__':
    main()

