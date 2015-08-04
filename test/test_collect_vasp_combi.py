#!/usr/bin/python
# -*- coding: utf-8 -*-
"""test"""
import os
import re
import glob
import pylab
import shutil
import unittest
import grapy
import collect_vasp
from commopy import Bash


def main():
    unittest.main()

TEST_PATH = '/Users/enoki/Researches/Analysis/Codes/01_testRun/'

class Collect(unittest.TestCase):
    path = os.path.join(TEST_PATH)

    def test_alt_elements_from_dir(self):
        path_list = glob.glob(os.path.join(self.path, 'combi_data', 'elem_*'))
        data = collect_vasp.Energy(path_list, 'POSCAR', 'OSZICAR')
        data.alt_elements_from_dir()

    def _test_correct_convex_hull_dat(self):
        path_list = glob.glob(os.path.join(self.path, 'Fe-Ni-Si', '*'))
        data = collect_vasp.Energy(path_list, 'POSCAR', 'OSZICAR')
        data.set_enthalpy()
        data.set_comp_dict_f()
        data.set_comp_dict_i()
        data.output_keys = ['comp_dict_f', 'enthalpy']
        print(data)

    def _test_correct_combi_enthalpy(self):
        path_list = glob.glob(os.path.join(self.path, 'combi_data', 'elem_*'))
        data = collect_vasp.Energy(path_list, 'POSCAR', 'OSZICAR')
        data.set_enthalpy()
        data.set_elements_z()  # 原子番号をset
        data.data.sort(key=lambda x: x['Z'][2])
        data.data.sort(key=lambda x: x['Z'][0])
        data.output_keys = ['elements', 'enthalpy']
        table_ene = data.table_combi([0, 2], 'enthalpy')
        print(table_ene)
        table_mag = data.table_combi([0, 2], 'mag')
        print(table_mag)
        table_cova = data.table_combi([0, 2], 'c/a')
        print(table_cova)

        table_ene = data.separate_data('elements', 0)
        print(table_ene[0]['elements'][:,2])
        plt = grapy.Vertical(3)
        self.plot(plt, table_ene, 'order', 'enthalpy', 'mag', 'c/a')

    def plot(self, plt, table, x, *y):
        for i in range(0, 3):
            plt.set_style(i*10+5, label=['test'])
            table[i].set_order()
            plt.set123(table[i], x, *y)
            pylab.xticks(table[i][x], table[i]['elements'][:,2])
        plt.ax1.legend(['B','C','N'], loc='lower right')
        plt.adjust_auto()
        pylab.xlim(-0.5, 20.5)
        plt.plot('show')

if __name__ == '__main__':
    main()
