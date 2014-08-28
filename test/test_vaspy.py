#!/usr/bin/python
# -*- coding: utf-8 -*-
"""test"""
import os
import glob
import vaspy
import shutil
import unittest
from commopy import Bash, Compare, Cabinet


def main():
    """main"""
    unittest.main()

TEST_PATH = ('/Users/enoki/Documents/01_ResearchData/Calculations/'
             '99_python/01_testRun/')

class VaspyOutcar(unittest.TestCase):
    """for OUTCAR"""
    path = os.path.join(TEST_PATH, 'vaspy', 'outcar')
    def test_outcar(self):
        """for outcar obj"""
        outcar = vaspy.Outcar(os.path.join(self.path, 'OUTCAR'))
        o_lines = Cabinet.read_file(os.path.join(self.path, 'OUTCAR'))
        print(outcar.get_mag(o_lines))

class VaspyIncar(unittest.TestCase): #pylint: disable=R0904
    """for Incar"""
    path = ('/Users/enoki/Documents/01_ResearchData/Calculations'
            '/99_python/01_testRun/prepVaspInputs/prepInputs/')
    vaspy.MakeInputs.all(path)
    read = vaspy.IncarReadPoscar.read_incar
    check_band = read(os.path.join(path, 'check_band'))
    check_soc001 = read(os.path.join(path, 'check_soc001'))
    check_volume = read(os.path.join(path, 'check_volume'))
    check_cell = read(os.path.join(path, 'check_cell'))

    out_band = read(os.path.join(path, 'INCAR_band'))
    out_soc001 = read(os.path.join(path, 'INCAR_soc001'))
    out_volume = read(os.path.join(path, 'INCAR_volume'))
    out_cell = read(os.path.join(path, 'INCAR_cell'))

    var1 = Compare.dict(check_band, out_band)
    var2 = Compare.dict(check_soc001, out_soc001)
    var3 = Compare.dict(check_volume, out_volume)
    var4 = Compare.dict(check_cell, out_cell)

    def test_makeall(self):
        """ for make all """
        print("VaspyIncar test_makeall")
        self.assertEqual(100, self.var1)
        self.assertEqual(100, self.var2)
        self.assertEqual(100, self.var3)
        self.assertEqual(100, self.var4)


class VaspyPoscar(unittest.TestCase): #pylint: disable=R0904
    """for Poscar"""
    path = ('/Users/enoki/Documents/01_ResearchData/Calculations'
            '/99_python/01_testRun/prepVaspInputs/series_vasp/')

    def test_template(self):
        """"template"""
        print("VaspyPoscar test_template")
        temp = vaspy.Poscar(None)
        self.assertEqual(temp.poscar_title, "FeCo (bench_test)\n")

    def test_readfrom_lines(self):
        """ read from lines format"""
        posc_lines = Cabinet.read_file(os.path.join(self.path, 'POSCAR'))
        poscar = vaspy.Poscar(posc_lines)
        print(poscar.poscar_title)

    def test_normalize(self):
        """normalize lattice parameter """
        print("VaspyPoscar test_normalize")
        os.chdir(self.path)
        clean_prev_dir('.', 'normalize')
        Bash.mkdir('normalize')
        poscar = vaspy.Poscar('POSCAR')
        poscar.normalize_lattice()
        poscar.write_poscar('normalize/POSCAR')
        poscar1 = vaspy.Poscar('POSCAR')
        poscar2 = vaspy.Poscar('normalize/POSCAR')
        len1 = poscar1.get_lattice_length()
        len2 = poscar2.get_lattice_length()
        vol1 = poscar1.get_cell_volume()
        vol2 = poscar1.get_cell_volume()
        self.assertEqual(len1, len2)
        self.assertEqual(vol1, vol2)
        self.assertEqual(poscar2.cell_lattices[0, 0], 1)
        self.assertEqual(poscar2.cell_lattices[0, 1], 0)
        self.assertEqual(poscar2.cell_lattices[0, 2], 0)


def clean_prev(path, files):
    """
    filesを消去
    """
    trush_list = glob.glob(os.path.join(path, files))
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
