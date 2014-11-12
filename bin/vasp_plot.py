#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
VASPのデータを収集してitxファイルを作成
"""
import os
import glob
import argparse

import itx
import collect_vasp

def main():
    """
    use argparse
    """
    parser = argparse.ArgumentParser(description='print series of data')
    # path (outputファイルまで入れて検索させる)
    parser.add_argument('strings', metavar='P', type=str, nargs='+',
                        help='enter path')
    # for FeB
    parser.add_argument('--feb', dest='run', const=for_uspex_FeB_print_data,
                        action='store_const', default=plot_mae_cova)
    #parser.add_argument('--def', dest='run', const=print_data)
    parser.add_argument('--fix_cova', dest='run', const=plot_mae_cova_02,
                        action='store_const', default=plot_mae_cova)
    args = parser.parse_args()
    print(args.strings)
    args.run(args.strings[0])


def plot_mae_cova(path):
    """
    energy, mag, meaのc/a依存性をプロット
    pathはoutputファイルまで検索させる
    """
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    data = collect_vasp.Energy(dir_list, 'POSCAR', 'OSZICAR')
    data.set_comp_dict_f()
    data.set_enthalpy()
    data.set_mae('OSZICAR_soc001', 'OSZICAR_soc100')
    data.data.sort(key=lambda x: x['c/a'])
    wavex = itx.Wave('cova', data['c/a'], r"\\Z20c/a")
    mu = r"\\F'Symbol'm\\F'Times New Roman'"
    wave0 = itx.Wave('mae', data['mae'], r"\\Z20MAE ({0}eV/atom)".format(mu))
    wave1 = itx.Wave('mag', data['mag'],
                     r"\\Z20Mag. ({0}\\M\\Z12B\\Z20/atom)".format(mu))
    wave2 = itx.Wave('enthalpy', data['enthalpy'], r"\\Z20Enthalpy (kJ/atom)")
    waves = [wavex, wave0, wave1, wave2]

    title = os.path.basename(os.path.abspath('.'))
    pref = itx.Produce(title, waves)
    pref.vertical3(waves)
    search_path = "/home/enoki/Dropbox/plot/{0}_0*.itx".format(title)
    num_list = [int(x.split('_')[-1].split('.')[0])
                for x in glob.glob(search_path)]
    num = 1
    if num_list:
        num = max(num_list) + 1
    fname = "/home/enoki/Dropbox/plot/{0}_{1:0>3}.itx".format(title, num)
    with open(fname, 'w') as wfile:
        wfile.write(pref.to_itx())

def plot_mae_cova_02(path):
    """
    energy, mag, meaのc/a依存性をプロット
    pathはoutputファイルまで検索させる
    primitive cellで軸の取り方がtetragonalでないケース
    ファイル名からc/aを取り直す
    """
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    data = collect_vasp.Energy(dir_list, 'POSCAR', 'OSZICAR')
    data.set_comp_dict_f()
    data.set_enthalpy()
    data.set_mae('OSZICAR_soc001', 'OSZICAR_soc100')
    # fix c/a
    for dat in data.data:
        dat['c/a'] = float(dat['path'].split('_')[-1])
    data.data.sort(key=lambda x: x['c/a'])
    wavex = itx.Wave('cova', data['c/a'], r"\\Z20c/a")
    mu = r"\\F'Symbol'm\\F'Times New Roman'"
    wave0 = itx.Wave('mae', data['mae'], r"\\Z20MAE ({0}eV/atom)".format(mu))
    wave1 = itx.Wave('mag', data['mag'],
                     r"\\Z20Mag. ({0}\\M\\Z12B\\Z20/atom)".format(mu))
    wave2 = itx.Wave('enthalpy', data['enthalpy'], r"\\Z20Enthalpy (kJ/atom)")
    waves = [wavex, wave0, wave1, wave2]

    title = os.path.basename(os.path.abspath('.'))
    pref = itx.Produce(title, waves)
    pref.vertical3(waves)
    search_path = "/home/enoki/Dropbox/plot/{0}_0*.itx".format(title)
    num_list = [int(x.split('_')[-1].split('.')[0])
                for x in glob.glob(search_path)]
    num = 1
    if num_list:
        num = max(num_list) + 1
    fname = "/home/enoki/Dropbox/plot/{0}_{1:0>3}.itx".format(title, num)
    with open(fname, 'w') as wfile:
        wfile.write(pref.to_itx())

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
