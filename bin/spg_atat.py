#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ATATのstr.outの空間群を求める
"""
import os
import numpy as np
from collections import defaultdict
import pymatgen as mg
from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifWriter
from pymatgen.symmetry.finder import SymmetryFinder

with open('str.out', 'r') as rfile:
    lines = rfile.read()
strs = lines.split('\nend\n\n')

def get_symetry(lines, output):
    latt_tmp = [[float(x) for x in y.split()] for y in lines[0:3]]
    trans = [[float(x) for x in y.split()] for y in lines[3:6]]
    sites = [[float(x) for x in y.split()[0:3]] for y in lines[6:]]
    elements = [x.split()[3] for x in lines[6:]]
    latt = np.dot(np.array(trans),np.array(latt_tmp))
    latt_obj = mg.Lattice(latt)
    structure = mg.Structure(latt, elements, sites)
    finder = SymmetryFinder(structure)
    out = finder.get_spacegroup_symbol() + "\n"
    #structure = finder.find_primitive()
    #structure = finder.get_primitive_standard_structure()
    poscar = Poscar(structure)
    #cif = CifWriter(structure)
    out += poscar.__str__()
    #print(cif)
    return out

# out = []
# for i in range(0, len(strs[:-1])):
#     out.append(get_symetry(strs[i].split('\n'), "strs/POSCAR.{0}".format(i)))
# print(len(out))
# setl = list(set(out))
# print(len(setl))
# for item in setl:
    # print(item)

def eight_sub_latt(structure):
    """
    8 sub latticeに属する構造の場合 Trueをretrun
    """
    elements = [x.split()[3] for x in structure.split('\n')[6:]]
    # if elements.count('Fe') % 2 != 0:
    #     return False
    # if elements.count('Fe') < 4:
    #     return False
    # if elements[0] != elements[7]:
    #     return False
    # if elements[1] != elements[6]:
    #     return False
    # if elements[2] != elements[5]:
    #     return False
    # if elements[3] != elements[4]:
    #     return False
    # if elements[8] != elements[15]:
    #     return False
    # if elements[9] != elements[14]:
    #     return False
    # if elements[10] != elements[13]:
    #     return False
    # if elements[11] != elements[12]:
    #     return False
    return True

judge = []
for i in range(0, len(strs[:-1])):
    judge.append(eight_sub_latt(strs[i]))


i = 0
out_str = []
for lines, boo in zip(strs[:-1], judge):
    if boo:
        i += 1
        lines = lines.split('\n')
        latt_tmp = [[float(x) for x in y.split()] for y in lines[0:3]]
        trans = [[float(x) for x in y.split()] for y in lines[3:6]]
        sites = [[float(x) for x in y.split()[0:3]] for y in lines[6:]]
        elements = [x.split()[3] for x in lines[6:]]
        latt = np.dot(np.array(trans), np.array(latt_tmp))
        latt_obj = mg.Lattice(latt)
        structure = mg.Structure(latt, elements, sites)
        finder = SymmetryFinder(structure)
        print(i)
        print(finder.get_spacegroup_symbol())
        structure = finder.get_primitive_standard_structure()
        out_str.append(Poscar(structure).__str__())
        print((Poscar(structure)))

#final_str = list(set(out_str))
final_str = out_str
print(len(final_str))

i = 0
for structure in final_str:
    i += 1
    posc = structure.split('\n')
    site_line = [" ".join(x.split()[0:3]) for x in posc[8:]]
    posc[8:] = site_line
    os.makedirs('str{0}'.format(i))
    with open('str{0}/POSCAR'.format(i), 'w') as wfile:
        for line in posc:
            wfile.write(line + '\n')
    # print(i)
    # for line in posc:
    #     print(line)
