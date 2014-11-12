#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""pymatgenのwrapper"""
import os
import pickle
from pymatgen.matproj.rest import MPRester
from pymatgen.io import vaspio

class Incar(vaspio.Incar):
    """
    ToDo outputの形式を見やすくする
    """
    pass

class Ternary(object):
    """
    三元系のデータを取り扱う
    """
    source_path = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
                   "99_python/00_workSpace/sorce/pymatgen/ternary_convex")
    def __init__(self, elements):
        mpr = MPRester("WTxsDhRV7g2Mcbqw")
        self.elements = elements
        self.entries = mpr.get_entries_in_chemsys(elements)

    def to_convex_hull(self, unit='kJ/mol'):
        """
        convex_hull用の形式でデータをreturn
        return: end_memb, not_end
        """
        factor = {'kJ/mol': 96.485344520851, 'eV/atom': 1}[unit]
        data_list = []
        for entry in self.entries:
            formula = entry.as_dict()['data']['unit_cell_formula']
            formation_e = entry.as_dict()['data']['formation_energy_per_atom']
            if formation_e <= 0:
                single_data = []
                sum_atoms = sum([x for x in formula.values()])
                for element in self.elements:
                    try:
                        single_data.append(formula[element]/sum_atoms)
                    except KeyError:
                        single_data.append(0)
                single_data.append(formation_e*factor)
                single_data.append(entry.as_dict()['data']['volume']/sum_atoms)
                data_list.append(single_data)
        data_list.sort(key=lambda x: x[3], reverse=True)
        end_memb = [[x[0], x[1], x[3], x[4]] for x in data_list[0:3]]
        not_end = [[x[0], x[1], x[3], x[4]] for x in data_list[3:]]
        fname = "".join(self.elements) + ".pickle"
        pickle_path = os.path.join(self.source_path, fname)
        with open(pickle_path, 'wb') as wbfile:
            pickle.dump([end_memb, not_end], wbfile)
        return end_memb, not_end





