#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Birch-Murnaghan Fitting
"""
import glob
import collect_vasp
from fitting_analysis import FitData


def main():
    """main"""
    #dir_list = glob.glob("volume_*")
    dir_list = glob.glob("QHA*/perfect")
    print_data(dir_list)

def print_data(dir_list):
    data = collect_vasp.Energy(dir_list, 'POSCAR', 'OSZICAR')
    data.data.sort(key=lambda x: x['volume'])
    data.output_keys = ['volume', 'energy']
    fit_res = FitData.Murnaghan_fit(data['volume']*8, data['energy']*8)
    for i in range(0, len(fit_res[0])):
        print("{0}\t{1}".format(fit_res[0][i], fit_res[1][i]))
    print(data)
    print(fit_res[2:])


if __name__ == '__main__':
    main()
