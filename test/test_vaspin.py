#!/usr/bin/python
# -*- coding: utf-8 -*-
"""test"""
import os
import glob
import vaspy
import shutil
import unittest
from commopy import Bash, Compare, Cabinet
from pymatgen.io import vaspio
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vaspio_set import MITVaspInputSet
import vaspin


def main():
    """main"""
    unittest.main()

TEST_PATH = ('/Users/enoki/Documents/01_ResearchData/Calculations/'
             '99_python/01_testRun/')

class TestVaspio(unittest.TestCase):
    """pymatgenの動作をチェック"""
    path = os.path.join(TEST_PATH, 'vaspin', 'pymatgen')
    def _testPotcar(self):
        """
        Potcarの動作チェック
        """
        potcar = vaspio.Potcar(symbols=['Fe', 'N', 'Co'])

        enmax = max([x.enmax for x in potcar])
        self.assertTrue(enmax == 400)

        rwigs = [x.rwigs for x in potcar]
        self.assertTrue(rwigs == [1.302, 0.741, 1.302])

        dst = os.path.join(self.path, 'POTCAR')
        potcar.write_file(dst)

    def _testPoscar(self):
        """
        Poscarをチェック
        """
        src = os.path.join(self.path, 'src_poscar')
        poscar = vaspio.Poscar.from_file(src, check_for_POTCAR=False)

        symbols = poscar.site_symbols
        self.assertTrue(symbols == ['Fe', 'Co', 'Ni'])
        print(dir(poscar.structure))
        print(poscar.structure.reciprocal_lattice.abc)
        new_symbols = ['Mn', 'Cr', 'Fe']
        print(poscar.structure.symbol_set)

    def _testKpoints(self):
        """
        Kpointsをチェック
        """
        # for Band
        src = os.path.join(self.path, 'src_poscar')
        poscar = vaspio.Poscar.from_file(src, check_for_POTCAR=False)
        hsk = HighSymmKpath(poscar.structure)
        kpts = hsk.get_kpoints()
        args = {'comment': "Kpoints for band calc",
                'kpts': kpts[0],
                'num_kpts': len(kpts[0]),
                'labels': kpts[1],
                'style': 'Reciprocal',
                'kpts_weights': [1]*len(kpts[0])}
        kpoints = vaspio.Kpoints(**args)
        print(kpoints)

    def testVaspio_set(self):
        #SOURCE_DIR = vasp_input_set

        src = os.path.join(self.path, 'src_poscar')
        dst = os.path.join(self.path, 'vaspio_set')
        poscar = vaspio.Poscar.from_file(src, check_for_POTCAR=False)
        #MITVaspInputSet().write_input(poscar.structure, dst)

        print(dir(MITVaspInputSet()))
        #print(MITVaspInputSet().potcar_settings)
        m = MITVaspInputSet()
        m.incar_settings['NSW'] = 10
        m.incar_settings['NELM'] = 60

        print(m.incar_settings)
        m.write_input(poscar.structure, dst)

class TestMyVaspin(unittest.TestCase):
    PATH = os.path.join(TEST_PATH, 'vaspin', 'myinputs', 'new')
    def test_relax(self):
        src = os.path.join(self.PATH, 'POSCAR_src')
        inputs = vaspin.ProduceInputs.from_file(src)
        os.chdir(self.PATH)
        os.chdir('relax_all')
        inputs.write_relax('all')
        os.chdir('../relax_volume')
        inputs.write_relax('volume')
        os.chdir('../relax_ion')
        inputs.write_relax('ion')

    def test_static(self):
        src = os.path.join(self.PATH, 'POSCAR_src')
        inputs = vaspin.ProduceInputs.from_file(src)
        os.chdir(self.PATH)
        os.chdir('static')
        inputs.write_static()

    def test_non_colli(self):
        src = os.path.join(self.PATH, 'POSCAR_src')
        inputs = vaspin.ProduceInputs.from_file(src)
        os.chdir(self.PATH)
        os.chdir('non_colli')
        inputs.write_non_collinear()











if __name__ == '__main__':
    main()
