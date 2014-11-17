#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
POSCARをpymatgenで取り扱う
"""
import os
import numpy as np
import glob
from collections import defaultdict
import pymatgen as mg
from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifWriter
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.analysis.structure_matcher import StructureMatcher

def get_primitive(fname):
    poscar = Poscar.from_file(fname)

    finder = SymmetryFinder(poscar.structure)
    spg_num = finder.get_spacegroup_number()

    primitive = finder.get_primitive_standard_structure()
    return spg_num, primitive

    # prim_pos = Poscar(primitive)
    # prim_pos.comment += "    (" + str(spg_num) + ": " + spg + ")"
    # prim_pos.write_file('POSCAR.prim')

    # finder = SymmetryFinder(primitive)
    # standard = finder.get_conventional_standard_structure()
    # stand_pos = Poscar(standard)
    # stand_pos.comment += "    (" + str(spg_num) + ": " + spg + ")"
    # stand_pos.write_file('POSCAR.stand')

    # match = StructureMatcher()
    # print(match.fit(standard, primitive))

def for_zengen(path):
    path_list = glob.glob(path)
    structures = defaultdict(list)
    for path in path_list:
        spg_num, primitive = get_primitive(path)
        structures[spg_num].append(primitive)
    irreps = []
    for spg_num in sorted(structures.keys(), reverse=True):
        while len(structures[spg_num]) != 1:
            struct0 = structures[spg_num].pop()
            matcher = StructureMatcher()
            judge = [matcher.fit(struct0, x) for x in structures[spg_num]]
            if not True in judge:
                irreps.append(struct0)
        irreps.append(structures[spg_num][0])
    return irreps

def produce(irreps):
    os.makedirs('irrep')
    for i, irrep in enumerate(irreps):
        poscar = Poscar(irrep)
        symbols = poscar.site_symbols
        natoms = poscar.natoms
        name_dict = {'Al': 'A', 'Ti': 'B'}
        tmp = ["{0}{1}".format(name_dict[x], y)
               for x, y in zip(symbols, natoms)]
        finder = SymmetryFinder(irrep)
        spg_num = finder.get_spacegroup_number()
        spg = "_".join(finder.get_spacegroup_symbol().split('/'))
        dirname = "No." + "{0:03d}".format(i) + "_" + spg + "_" + "".join(tmp)
        poscar.comment += "    (#" + str(spg_num) + ": " + spg + ")"

        standard = finder.get_conventional_standard_structure()
        stand_pos = Poscar(standard)
        stand_pos.comment += "    (#" + str(spg_num) + ": " + spg + ")"

        os.makedirs(os.path.join('irrep', dirname))
        poscar.write_file(os.path.join('irrep', dirname, 'POSCAR.prim'))
        stand_pos.write_file(os.path.join('irrep', dirname, 'POSCAR.std'))



if __name__ == '__main__':
    irreps = for_zengen('*/POSCAR.ini')
    produce(irreps)
