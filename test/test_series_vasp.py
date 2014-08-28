#!/usr/bin/python
# -*- coding: utf-8 -*-
"""test"""
import os
import copy
import random
import shutil
import unittest
import series_vasp
import vaspy
from test_judgeRX import VaspyIncar
from commopy import Bash, Compare


def main():
#    a = test_judgeRX.VaspyIncar()
    unittest.main()


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
    trush_list = Bash.find_files(path, dirc)
    for trush in trush_list:
        fname = os.path.join(path, trush)
        shutil.rmtree(fname)
        print("{0} is removed.".format(fname))


class SeriesVasp(unittest.TestCase):
    path = ('/Users/enoki/Documents/01_ResearchData/Calculations'
            '/99_python/01_testRun/prepVaspInputs/series_vasp/')

    def test_cova(self):
        print("SeriesVasp test_cova")
        os.chdir(self.path)
        clean_prev_dir('.', 'cova')
        Bash.mkdir('cova')
        cova_list = [round(random.uniform(1, 3), 3) for x in range(0, 3)]
        series = series_vasp.Produce('POSCAR', 'cova')
        series.set_cova(cova_list)
        series.make_files()

        orig_pos = vaspy.Poscar('POSCAR')
        vol_orig = orig_pos.get_cell_volume()

        test_param = random.choice(series.series)
        test_path = test_param['path']
        test = vaspy.Poscar(os.path.join(test_path, 'POSCAR'))

        len_a, _, len_c = test.get_lattice_length()
        vol = test.get_cell_volume()
        self.assertAlmostEqual(len_c/len_a, test_param['cova'])
        self.assertAlmostEqual(vol_orig, vol)

    def test_volume(self):
        print("SeriesVasp test_volume")
        os.chdir(self.path)
        clean_prev_dir('.', 'volume')
        Bash.mkdir('volume')
        volume_list = [5, 6, 7]
        series = series_vasp.Produce('POSCAR', 'volume')
        series.set_volume(volume_list)
        series.make_files()

    def test_latt(self):
        print("SeriesVasp test_latt")
        os.chdir(self.path)
        clean_prev_dir('.', 'scale')
        Bash.mkdir('scale')
        scale_list = [round(random.uniform(1, 3), 3) for x in range(0, 3)]
        series = series_vasp.Produce('POSCAR', 'scale')
        series.set_scale(scale_list)
        series.make_files()

        test_param = random.choice(series.series)
        test_path = test_param['path']
        test = vaspy.Poscar(os.path.join(test_path, 'POSCAR'))

        orig_pos = vaspy.Poscar('POSCAR')
        orig_a, _, orig_c = orig_pos.get_lattice_length()

        len_a, _, len_c = test.get_lattice_length()
        scale = test.cell_scale
        self.assertAlmostEqual(len_c/len_a, orig_c/orig_a)
        self.assertAlmostEqual(scale, test_param['scale'])

    def test_double(self):
        print("SeriesVasp test_double")
        os.chdir(self.path)
        clean_prev_dir('.', 'scale_cova')
        scale_list = [round(random.uniform(1, 3), 3) for x in range(0, 3)]
        cova_list = [round(random.uniform(1, 3), 3) for x in range(0, 3)]

        series = series_vasp.Produce('POSCAR', 'scale_cova')
        series.set_cova(cova_list)
        series.set_scale(scale_list)

        series.make_files()

        test_param = random.choice(series.series)
        test_path = test_param['path']
        test = vaspy.Poscar(os.path.join(test_path, 'POSCAR'))
        len_a, _, len_c = test.get_lattice_length()
        scale = test.cell_scale
        self.assertAlmostEqual(len_c/len_a, test_param['cova'])
        self.assertAlmostEqual(scale, test_param['scale'])

    def test_incar_tag(self):
        print("SeriesVasp test_incar_tag")
        os.chdir(self.path)
        clean_prev_dir('.', 'incar_tag')
        series = series_vasp.Produce('POSCAR', 'incar_tag')
        tag_list = [100, 200]
        series.set_incar_tag('encut', tag_list)
        series.make_files()
        path = series.series[0]['path']
        read = vaspy.IncarReadPoscar.read_incar
        check_100 = read(os.path.join(path, 'INCAR_cell'))
        self.assertEqual(check_100['encut'], 130)

    def test_incar_fixedtag(self):
        print("SeriesVasp test_incar_fixedtag")
        os.chdir(self.path)
        clean_prev_dir('.', 'incar_fixedtag')
        series = series_vasp.Produce('POSCAR', 'incar_fixedtag')
        tag_list = [100, 200]
        series.set_incar_fixedtag('encut', tag_list)
        series.make_files()
        series.make_list_run("run.sh")
        series.append_list_run("run.sh")

        path = series.series[0]['path']
        read = vaspy.IncarReadPoscar.read_incar
        check_100 = read(os.path.join(path, 'INCAR_cell'))
        self.assertEqual(check_100['encut'], 100)

    def test_incar_altmagmom(self):
        print("SeriesVasp incar_altmagmom")
        os.chdir(self.path)
        clean_prev_dir('.', 'incar_altmagmom')
        vaspy.IncarReadPoscar.cls_add_extratag({'magmom': [2, -2, 1]})
        series = series_vasp.Produce('POSCAR', 'incar_altmagmom')
        tag_list = [100, 200]
        series.set_incar_tag('encut', tag_list)
        series.make_files()

        path = series.series[0]['path']
        read = vaspy.IncarReadPoscar.read_incar
        check = read(os.path.join(path, 'INCAR_cell'))
        magmom = check['magmom']
        self.assertTrue(magmom == [2, -2, 1])
        vaspy.IncarReadPoscar.cls_initialize()

    def test_elements(self):
        print("SeriesVasp elements")
        os.chdir(self.path)
        clean_prev_dir('.', 'elements')
        series = series_vasp.Produce('POSCAR', 'elements')
        elements_list = [['Co', 'Sb'], ['Ni', 'Si'], ['Pb', 'Al']]
        series.set_elements(elements_list)
        series.make_files()

if __name__ == '__main__':
    main()
