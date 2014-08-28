#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
mae計算
"""
import collect_vasp
import glob
import os

def main():
    """main"""
    dir_list = glob.glob("EA*/OSZICAR")
    dir_list = [os.path.dirname(x) for x in dir_list]
    print_data(dir_list)


def print_data(dir_list):
    os.chdir('results1_POSCARS_FeB')
    dir_list = glob.glob("EA*/OSZICAR")
    dir_list = [os.path.dirname(x) for x in dir_list]
    data = collect_vasp.Energy(dir_list, 'POSCAR', 'OSZICAR')
    data.set_enthalpy()
    data.set_comp_dict_f()
    data.set_comp_dict_i()
    data.output_keys = ['comp_dict_f', 'path', 'enthalpy']
    new = {}
    new3 = {}
    for data in data.data:
        new.update({data['path']: data['enthalpy']})
        new3.update({data['path']: data['comp_dict_i']})

    os.chdir('../results1_POSCARS_Fe_only')
    dir_list = glob.glob("EA*/OSZICAR")
    dir_list = [os.path.dirname(x) for x in dir_list]
    data = collect_vasp.Energy(dir_list, 'POSCAR', 'OSZICAR')
    data.set_enthalpy()
    data.set_comp_dict_f()
    data.set_comp_dict_i()
    data.output_keys = ['comp_dict_f', 'path', 'enthalpy']
    new2 = {}
    for data in data.data:
        new2.update({data['path']: data['enthalpy']})
    lines = ""
    for key in new:
        try:
            e = new[key] - new2[key]
            lines += " ".join([key, str(e), str(new3[key]['B'])]) + "\n"
        except KeyError:
            pass
    print(lines)


if __name__ == '__main__':
    main()
