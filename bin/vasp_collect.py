#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
VASPのデータを収集してprint
"""
import os
import glob
import argparse
import collect_vasp

def main():
    """
    use argparse
    """
    parser = argparse.ArgumentParser(description='print series of data')
    parser.add_argument('strings', metavar='P', type=str, nargs='+',
                        help='enter path')
    parser.add_argument('--feb', dest='run', const=for_uspex_FeB_print_data,
                        action='store_const', default=print_data)
    #parser.add_argument('--def', dest='run', const=print_data)
    args = parser.parse_args()
    print(args.strings)
    args.run(args.strings[0])


def print_data(path):
    """
    入力のpathはワイルドカードを使用可能
    os.path.dirnameを修得するのでoutputファイルまで入れて検索させること
    """
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    data = collect_vasp.Energy(dir_list, 'POSCAR', 'OSZICAR')
    #data.set_enthalpy()
    print(data)

def for_uspex_FeB_print_data():
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
