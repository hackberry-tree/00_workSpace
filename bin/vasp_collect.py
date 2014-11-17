#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
VASPのデータを収集してprint
"""
import os
import glob
import math
import pylab
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
    parser.add_argument('strings', metavar='P', type=str, nargs='+',
                        help='enter path')
    # for FeB
    parser.add_argument('--feb', dest='run', const=for_uspex_FeB_print_data,
                        action='store_const', default=print_data)
    #parser.add_argument('--def', dest='run', const=print_data)
    parser.add_argument('--mae', dest='mae', const=True, action='store_const',
                        default=False)
    args = parser.parse_args()
    print(args.strings)
    args.run(args.strings[0], args.mae)


def print_data(path, mae='False'):
    """
    入力のpathはワイルドカードを使用可能
    os.path.dirnameを修得するのでoutputファイルまで入れて検索させること
    """
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    data = collect_vasp.Energy(dir_list, 'POSCAR', 'OSZICAR')
    data.set_comp_dict_f()
    data.set_enthalpy()
    if mae:
        data.set_mae('OSZICAR_soc001', 'OSZICAR_soc100')
    #data.output_keys = ['comp_dict_f', 'energy', 'enthalpy', 'path']
    print(data)

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
            out[dirname].update({'comp_dict_f': data['comp_dict_f']})
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
