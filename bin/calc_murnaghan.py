#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Birch-Murnaghan Fitting
"""
import os
import sys
import glob
import collect_vasp
import numpy as np
import pylab
from fitting_analysis import FitData
from pymatgen.io.vaspio import Poscar

def main():
    """main"""
    if len(sys.argv) >= 2:
        get_results(sys.argv[1])
        return
    get_result_single('.')

def get_results(path):
    """
    入力のpathはワイルドカードを使用可能
    os.path.dirnameを修得するのでoutputファイルまで入れて検索させること
    e.g. "GS*/finish_voldep"
    """
    dir_list = [os.path.dirname(x) for x in glob.glob(path)]
    for dirc in dir_list:
        # get_result_single(dirc)
        get_result_with_mag(dirc)


def get_result_single(path):
    """
    Murnaghan で fit した eps と 結果のファイルを作成する
    """
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
    fig = pylab.figure()
    pylab.plot(data['volume'], data['energy'])
    pylab.plot(fit_res[0], fit_res[1])

    judge = judge_fit(fit_res, data)
    with open(os.path.join(path, 'is_fit_ok'), 'w') as wfile:
        wfile.write(judge)

    lines = ""
    lines += "E:{0[2][0]} dE:{0[3][0]}\n".format(fit_res)
    lines += "B:{0[2][1]} dB:{0[3][1]}\n".format(fit_res)
    lines += "B':{0[2][2]} dB':{0[3][2]}\n".format(fit_res)
    lines += "V:{0[2][3]} dV:{0[3][3]}".format(fit_res)
    fig.text(0.15, 0.75, lines)
    fname = os.path.join(path, 'murnaghan_plot.eps')
    fig.savefig(fname)
    pylab.close('all')

    lines = ""
    lines += "Energy\t{0[2][0]}\terrE\t{0[3][0]}\n".format(fit_res)
    lines += "B\t{0[2][1]}\terrB\t{0[3][1]}\n".format(fit_res)
    lines += "B1\t{0[2][2]}\terrB1\t{0[3][2]}\n".format(fit_res)
    lines += "Volume\t{0[2][3]}\terrV\t{0[3][3]}\n".format(fit_res)
    fname = os.path.join(path, 'fit_results.dat')
    with open(fname, 'w') as wfile:
        wfile.write(lines)
    return fit_res


def get_result_with_mag(path):
    """
    Murnaghan で fit した eps と 結果のファイルを作成する
    total の磁気モーメントもプロットする
    """
    osz_path = os.path.join(path, 'voldep/volume_*/OSZICAR')
    dir_list = [os.path.dirname(x) for x in glob.glob(osz_path)]
    data = collect_vasp.Energy(dir_list, 'POSCAR', 'OSZICAR')
    data.data.sort(key=lambda x: x['volume'])
    data.set_mag_each_site()
    data.output_keys = ['symbol_mag/site', 'mag/site']
    mag = np.array(data['mag/site'])
    fig = pylab.figure(figsize=(6, 9))
    ax = pylab.subplot(2, 1, 2)

    for i in range(len(data.data[0]['symbol_mag/site'])):
        ax.plot(data['volume'], mag[:, i],
                   label=data.data[0]['symbol_mag/site'][i])
    ax.legend(loc='upper left', prop={'size':8})

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
    ax2 = pylab.subplot(2, 1, 1)
    ax2.plot(data['volume'], data['energy'])
    ax2.plot(fit_res[0], fit_res[1])

    judge = judge_fit(fit_res, data)
    with open(os.path.join(path, 'is_fit_ok'), 'w') as wfile:
        wfile.write(judge)

    lines = ""
    lines += "E:{0[2][0]} dE:{0[3][0]}\n".format(fit_res)
    lines += "B:{0[2][1]} dB:{0[3][1]}\n".format(fit_res)
    lines += "B':{0[2][2]} dB':{0[3][2]}\n".format(fit_res)
    lines += "V:{0[2][3]} dV:{0[3][3]}".format(fit_res)
    fig.text(0.15, 0.8, lines)
    fname = os.path.join(path, 'murnaghan_plot.eps')
    fig.savefig(fname)
    pylab.close('all')

    lines = ""
    lines += "Energy\t{0[2][0]}\terrE\t{0[3][0]}\n".format(fit_res)
    lines += "B\t{0[2][1]}\terrB\t{0[3][1]}\n".format(fit_res)
    lines += "B1\t{0[2][2]}\terrB1\t{0[3][2]}\n".format(fit_res)
    lines += "Volume\t{0[2][3]}\terrV\t{0[3][3]}\n".format(fit_res)
    fname = os.path.join(path, 'fit_results.dat')
    with open(fname, 'w') as wfile:
        wfile.write(lines)
    return fit_res


def judge_fit(fit_res, data):
    """
    range が over している結果
    誤差が大きい結果を拾う
    また誤差に 0 が与えられている結果も error とする
    誤差が大きいものは fitting そのものがうまくいっていない可能性があるので
    range よりも優先して出力する (上書きする)
    """
    judge = ""

    if fit_res[2][3] < data.data[0]['volume']:
        print('volume-range is too_large')
        judge = "too-large-range\n"

    if fit_res[2][3] > data.data[-1]['volume']:
        print('volume-range is too_small')
        judge = "too-small-range\n"

    if fit_res[3][0] == 0:
        print('error is too large')
        judge = "too-large-error\n"

    if fit_res[3][0] > 0.05:
        print('error is too large')
        judge = "too-large-error\n"

    if fit_res[3][2] > 3.00:
        print('error is too large')
        judge = "too-large-error\n"

    # if fit_res[3][3] > 0.02:
    #     print('error is too large')
    #     judge = "too-large-error\n"

    return judge


if __name__ == '__main__':
    main()
