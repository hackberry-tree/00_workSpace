#!/opt/local/bin/python
# -*- coding: utf-8 -*-
"""
energy, c/a, mag, volumeなどの情報を読んでdata_arrayを作成する
"""
import os
import glob
import pylab
import shutil
import unittest
import grapy
import collect_vasp
from commopy import Bash, Cabinet
from makeSeries import Combinatorial


def main():
    unittest.main()

TEST_PATH = ('/Users/enoki/Documents/01_ResearchData/Calculations/'
             '99_python/01_testRun/')


class Grapy(unittest.TestCase):
    path = os.path.join(TEST_PATH, 'grapy_Heusler')

    def test_plt(self):
        clean_prev_dir(self.path, 'output')
        Bash.mkdir(os.path.join(self.path, 'output', 'text'))
        combi_para = [['Co'], ['Al']]
        pathA = os.path.join(self.path, 'typeA')
        pathB = os.path.join(self.path, 'typeB')
        path_out = os.path.join(self.path, 'output')
        self.plt_double(combi_para, pathA, pathB, path_out)

    def plt_double(self, combipara, in_dir1, in_dir2, out_dir):
        """2種の結果のlattice依存性をプロット"""
        combi = Combinatorial(*combipara)
        combi.set_formula(2, 1)
        for composition in combi.compositions:
            formula = composition['formula']
            formula = ''.join(formula.split('2'))
            total_elem = ['Fe'] + composition['elements']
            total_num_atoms = [1] + composition['num_atoms']
            print(formula)

            all_dir_list = glob.glob(os.path.join(in_dir1, formula, 'fixed_*'))
            regular = collect_vasp.Energy(all_dir_list, 'POSCAR', 'OSZICAR')
            regular.data['c/a'] /= 2 ** 0.5
            regular.get_enthalpy(total_elem, total_num_atoms)
            regular.get_mae('OSZICAR_SOC001', 'OSZICAR_SOC100', 'mae')

            all_dir_list = glob.glob(os.path.join(in_dir2, formula, 'fixed_*'))
            inverse = collect_vasp.Energy(all_dir_list, 'POSCAR', 'OSZICAR')
            inverse.data['c/a'] /= 2 ** 0.5
            inverse.get_enthalpy(total_elem, total_num_atoms)
            inverse.get_mae('OSZICAR_SOC001', 'OSZICAR_SOC100', 'mae')

            rlines = [str(regular)]
            ilines = [str(inverse)]
            fname_reg = "Fe2{0[0]}{0[1]}_reg.txt".format(composition['elements'])
            fname_inv = "Fe2{0[0]}{0[1]}_inv.txt".format(composition['elements'])
            Cabinet.write_file(os.path.join(out_dir, 'text', fname_reg), rlines)
            Cabinet.write_file(os.path.join(out_dir, 'text', fname_inv), ilines)


            def plot():
                plt = grapy.Vertical(3)
                plt.set_title("Heusler Fe$_2${0[0]}{0[1]}"
                                .format(composition['elements']))
                plt.set_style('blue')
                plt.set123(regular.data, 'c/a', 'enthalpy', 'mag', 'mae')
                plt.set_style('magenta')
                plt.set123(inverse.data, 'c/a', 'enthalpy', 'mag', 'mae')

                plt.adjust_auto()
                plt.ax2.set_ylim(-0.5, 3)
                plt.plot('show')
                #fname = "Fe2{0[0]}{0[1]}.eps".format(composition['elements'])
                #plt.plot(os.path.join(out_dir, fname))
            plot()


def clean_prev(path, files):
    """
    filesを消去
    """
    trush_list = Bash.find_files(path, files)
    for trush in trush_list:
        fname = os.path.join(path, trush)
        os.remove(fname)
        print("{0} is removed.".format(fname))


def clean_prev_dir(path, dirc):
    """
    dirctoryを消去
    """
    trush_list = glob.glob(os.path.join(path, dirc))
    for trush in trush_list:
        shutil.rmtree(trush)
        print("{0} is removed.".format(trush))


if __name__ == '__main__':
    main()
