#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Birch-Murnaghan Fitting
"""
import os
import glob
import collect_vasp
import pylab
from fitting_analysis import FitData


def main():
    """main"""
    dir_list = glob.glob("volume_*/OSZICAR")
    dir_list = [os.path.dirname(x) for x in dir_list]
    #dir_list = glob.glob("QHA*/perfect")
    get_result('.')

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
