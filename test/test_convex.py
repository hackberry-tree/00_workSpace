#!/usr/bin/python
# -*- coding: utf-8 -*-
"""test"""
import numpy as np
import os
import math
import pylab
import shutil
import unittest
import convex_hull
from convex_hull import FindGS
from commopy import Bash
from mpl_toolkits.mplot3d import Axes3D

from pymatgen.matproj.rest import MPRester
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter
from pymatgen.io.cifio import CifParser
from pymatgen.io.vaspio import Poscar


TEST_PATH = "/Users/enoki/Researches/Analysis/Codes/01_testRun/"

def main():
    unittest.main()


class TestConvex(unittest.TestCase):
    path = os.path.join(TEST_PATH, 'convex_hull', 'FeNiSi')

    def _test_triangle(self):
        point1 = [0, 1, 2]
        point2 = [1, 2, 3]
        point3 = [1, 1, 1]
        point4 = [2, 3, 4]
        point5 = [1, 1, 1]
        point6 = [2, 4, 6]

        triangleA = [point1, point2, point3]
        triangleB = [point4, point5, point6]
        triangleC = [point1, point5, point4]

        #FindGS.resplit_triangle(triangleA, triangleB)
        tri1, tri2, judge = FindGS.resplit_triangle(triangleA, triangleC)
        print(tri1)
        #FindGS.resplit_triangle(triangleA, triangleA)

        # tris = [triangleB, triangleC]
        # print(FindGS.find_adjacent_triangles(triangleA, tris))

    def test_select_base(self):
        """
        桁数を揃えておかないと無限ループに陥る
        """
        print("Test select bases")
        end_memb = [[0, 0, 0], [0, 1, 0], [1, 0, 0]]
        not_end = [[0.25, 0.75, -27.853794735939683], [0.25, 0.25, -41.16333894885379]]
        #not_end = [[0.3333333333, 0.3333333333, -11.6165778644]]
        bases = FindGS.collect_base_triangles(end_memb, not_end)


    def test_icvm(self):
        """
        icvm の test
        """
        path = os.path.join(TEST_PATH, 'convex_hull', 'icvm')
        data = convex_hull.iCVM_energies.from_file(
            os.path.join(path, 'energies.txt'))

        end1 = \
            (data.fract_per_atom[:, [0, 2, 3]] == [3/4, 1/4, 0]).prod(axis=-1) == 1
        end2 = \
            (data.fract_per_atom[:, [0, 2, 3]] == [3/4, 0, 1/4]).prod(axis=-1) == 1
        print(data.enthalpy[end1])
        print(data.enthalpy[end2])

        data_xyz = np.c_[data.fract_per_atom[:, 2], data.fract_per_atom[:, 3],
                         data.enthalpy]
        fig = pylab.figure()
        ax = Axes3D(fig)
        pt3d = convex_hull.PlotTriangularCoord(ax)
        pt3d.plt_dot(data_xyz)
        pylab.show()

        # no_int = (data.num_atoms[:, 0] == 0)
        # data_xy = np.c_[data.fract_per_atom[no_int, 2], data.enthalpy[[no_int]]]
        # print(data_xy)
        # pylab.plot(data_xy[:, 0], data_xy[:, 1], 'd')
        # pylab.show()



    def _test_draw_convex_hull(self):
        data = pylab.loadtxt(self.path, comments='#')
        initial_base = [list(x) for x in data[0:3, [0, 2, 3]]]
        not_base = [list(x) for x in data[3:, [0, 2, 3]]]

        fig = pylab.figure()
        ax = Axes3D(fig)
        convex_hull.draw_convex_hull(ax, initial_base, not_base,
                                     ['Fe', 'Ni', 'Si'], [-60, 5])
        #pint(not_base)
        pylab.show()

    def _test_from_pm(self):
        """
        MPのデータをプロット
        """
        ent_list = []
        mpr = MPRester("WTxsDhRV7g2Mcbqw")
        composition = ['Fe', 'Ni', 'Si']
        entries = mpr.get_entries_in_chemsys(composition)
        for entry in entries:
            formula = entry.as_dict()['data']['unit_cell_formula']
            formation_e = entry.as_dict()['data']['formation_energy_per_atom']
            if formation_e <= 0:
                single_data = []
                sum_atoms = 0
                for elements, index in formula.items():
                    sum_atoms += index
                for element in composition:
                    try:
                        single_data.append(formula[element] / sum_atoms)
                    except KeyError:
                        single_data.append(0)
                single_data.append(formation_e * 96.485344520851)
                ent_list.append(single_data)
        initial_base = [[x[0], x[1], x[3]] for x in ent_list[0:3]]
        not_base = [[x[0], x[1], x[3]] for x in ent_list[3:]]
        fig = pylab.figure()
        ax = Axes3D(fig)
        convex_hull.draw_convex_hull(ax, initial_base, not_base,
                                     ['Fe', 'Ni', 'Si'], [-60, 5],
                                     color='magenta')
        #pint(not_base)
        pylab.show()


def clean_prev(path, files):
    """
    既存filesを消去
    """
    trush_list = Bash.find_files(path, files)
    for trush in trush_list:
        fname = os.path.join(path, trush)
        os.remove(fname)
        print("{0} is removed.".format(fname))


def clean_prev_dir(path, dirc):
    """
    既存dirctoryを消去
    """
    trush_list = Bash.find_files(path, dirc)
    for trush in trush_list:
        fname = os.path.join(path, trush)
        shutil.rmtree(fname)
        print("{0} is removed.".format(fname))


if __name__ == '__main__':
    main()
