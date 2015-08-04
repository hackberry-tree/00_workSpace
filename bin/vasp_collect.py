#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
VASPのデータを収集してprint
"""
import os
import glob
import math
import pylab
import pickle
import argparse
import collect_vasp
from fitting_analysis import FitData
from commopy import DataBox

def main():
    """
    use argparse
    """
    parser = argparse.ArgumentParser(description='print series of data')
    # path (outputファイルまで入れて検索させる)
    parser.add_argument('path', type=str, nargs=1,
                        help='enter path')
    # for FeB
    parser.add_argument('--feb', dest='run', const=for_uspex_FeB_print_data,
                        action='store_const', default=print_data)
    # for ATAT
    parser.add_argument('--atat', dest='run', const=print_data_atat,
                        action='store_const', default=print_data)

    # for MAE
    parser.add_argument('--mae', dest='mae', const=True,
                        action='store_const', default=False)

    # for combi
    # e.g. "elem_*/OSZICAR" --combi --combi_pos 0 1
    parser.add_argument('--combi', dest='combi', type=int, nargs=2,
                        action='append', default=None)

    # for arb. file name
    parser.add_argument('--posc', dest='posc', type=str, nargs=1,
                        action='store', default=None)
    parser.add_argument('--osz', dest='osz', type=str, nargs=1,
                        action='store', default=None)

    parser.add_argument('--mg_fit', dest='run', const=print_data_mg,
                        action='store_const', default=print_data)



    args = parser.parse_args()
    opt = {'path': args.path[0]}
    if args.mae:
        opt.update({'mae': True})
    if args.combi:
        opt.update({'pos': args.combi[0]})
        args.run = correct_combi
    if args.posc:
        opt.update({'posc': args.posc[0]})
    if args.osz:
        opt.update({'osz': args.osz[0]})
    # print(opt)
    args.run(**opt)

def print_data_mg(path, mae=False):
    """
    MurnaghanFit の結果を収集する
    """
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    data = collect_vasp.Mg_fit_results.from_dir_list(dir_list)
    data.set_comp_dict_f()
    print(data)
    data.make_energies_txt()

def correct_combi(path, pos, mae=False):
    """
    入力のpathはワイルドカードを使用可能
    os.path.dirnameを修得するのでoutputファイルまで入れて検索させること
    e.g. "elem_*/OSZICAR"
    読み込むファイルは osz.static
    """
    print(pos)
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    data = collect_vasp.Energy(dir_list, 'POSCAR', 'OSZICAR')
    data.set_enthalpy()
    data.alt_elements_from_dir()
    data.alt_cova_posstsd()
    data.set_elements_z()  # 原子番号をset

    if mae:
        data.set_mae('OSZICAR_soc001', 'OSZICAR_soc100')

    print(data.data[0]['Z'])
    data.data.sort(key=lambda x: x['Z'][pos[0]])
    data.data.sort(key=lambda x: x['Z'][pos[1]])
    data.output_keys = ['elements', 'enthalpy']

    if mae:
        table_mae = data.table_combi(pos, 'mae')
        print(table_mae)

    table_ene = data.table_combi(pos, 'enthalpy')
    print(table_ene)
    table_mag = data.table_combi(pos, 'mag')
    print(table_mag)
    table_cova = data.table_combi(pos, 'c/a')
    print(table_cova)

    table_ene = data.separate_data('elements', pos[0])
    # print(table_ene[pos[0]]['elements'][:,pos[1]])
    #plt = grapy.Vertical(3)
    #self.plot(plt, table_ene, 'order', 'enthalpy', 'mag', 'c/a')


def print_data(path, mae=False, posc='POSCAR', osz='OSZICAR'):
    """
    入力のpathはワイルドカードを使用可能
    os.path.dirnameを修得するのでoutputファイルまで入れて検索させること
    """
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    data = collect_vasp.Energy(dir_list, posc, osz)
    data.set_comp_dict_f()
    data.set_enthalpy()
    print(mae)
    if mae == True:
        data.set_mae('OSZICAR_soc001', 'OSZICAR_soc100')
        data.output_keys = ['energy', 'enthalpy', 'path', 'mae']
    print(data)

def print_data_atat(path):
    """
    入力のpathはワイルドカードを使用可能
    os.path.dirnameを修得するのでoutputファイルまで入れて検索させること
    150615
    convex hull 用に data の pickle を保存するように修正
    """
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    data = collect_vasp.Energy(dir_list, 'POSCAR.static', 'OSZICAR.static')
    data.set_comp_dict_f()
    data.set_enthalpy()
    data.output_keys = ['comp_dict_f', 'energy', 'enthalpy', 'spg', 'spg_num',
                        'volume', 'mag', 'path']
    print(data)
    with open("collect_data.pickle", 'wb') as wbfile:
        pickle.dump(data, wbfile)

def for_uspex_FeB_print_data(_, _2):
    two_dir = ['results1_POSCARS_FeB', 'results1_POSCARS_Fe_only']
    #os.chdir('results1_POSCARS_FeB')
    def get_res(path):
        osz_path = os.path.join(path, 'EA*/OSZICAR')
        dir_list = [os.path.dirname(x) for x in glob.glob(osz_path)]
        data = collect_vasp.MurnaghanFit(dir_list)
        data.set_enthalpy()
        data.set_comp_dict_f()
        data.set_comp_dict_i()
        data.output_keys = ['comp_dict_f', 'path', 'enthalpy']
        out = {}
        for data in data.data:
            dirname = os.path.basename(data['path'])
            out.update({dirname: {'enthalpy': data['enthalpy']}})
            out[dirname].update({'comp_dict_f': data['comp_dict_f']}    )
            out[dirname].update({'comp_dict_i': data['comp_dict_i']})
            out[dirname].update({'errE': data['errE']*96.4})
        return out
    feb = get_res(two_dir[0])
    fe = get_res(two_dir[1])
    data_box = []
    for key in feb:
        try:
            e = feb[key]['enthalpy'] - fe[key]['enthalpy']
            err = ((feb[key]['errE']**2) + (fe[key]['errE']**2))**0.5
            e2 = feb[key]['enthalpy'] * \
                 sum(feb[key]['comp_dict_i'].values())
            e2 -= fe[key]['enthalpy'] * \
                  sum(fe[key]['comp_dict_i'].values())
            data_box.append({'path': key, 'enthalpy': e, 'total_e': e2,
                                'comp_dict_f': feb[key]['comp_dict_f'],
                                'comp_dict_i': feb[key]['comp_dict_i'],
                                'errE': err})
        except KeyError:
            pass
    result = DataBox(data_box)
    print(result)

def get_result(path):
    osz_path = os.path.join(path, 'voldep/volume_*/OSZICAR')
    dir_list = [os.path.dirname(x) for x in glob.glob(osz_path)]
    data = collect_vasp.Energy(dir_list, 'POSCAR', 'OSZICAR')
    data.data.sort(key=lambda x: x['volume'])
    data.output_keys = ['volume', 'energy']
    fit_res = FitData.Murnaghan_fit(data['volume'], data['energy'])
    lines = ""
    for i in range(0, len(fit_res[0])):
        lines += "{0}\t{1}\n".format(fit_res[0][i], fit_res[1][i])
    fname = os.path.join(path, 'fit_curve.dat')
    with open(fname, 'w') as wfile:
        wfile.write(lines)
    pylab.plot(data['volume'], data['energy'])
    pylab.plot(fit_res[0], fit_res[1])
    fname = os.path.join(path, 'murnaghan_plot.eps')
    pylab.savefig(fname)
    lines = ""
    lines += "Energy\t{0[2][0]}\terrE\t{0[3][0]}\n".format(fit_res)
    lines += "B\t{0[2][1]}\terrB\t{0[3][1]}\n".format(fit_res)
    lines += "B1\t{0[2][2]}\terrB1\t{0[3][2]}\n".format(fit_res)
    lines += "Volume\t{0[2][3]}\terrV\t{0[3][3]}\n".format(fit_res)
    fname = os.path.join(path, 'fit_results.dat')
    with open(fname, 'w') as wfile:
        wfile.write(lines)
    return fit_res


if __name__ == '__main__':
    main()
