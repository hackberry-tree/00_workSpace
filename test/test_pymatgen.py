#!/usr/bin/env python
import pymatgen as mg

# random構造の作成
def random():
    lattice = mg.Lattice.cubic(2.78)
    structure = mg.Structure(lattice, ['Al'], [[0, 0, 0]])
    structure.make_supercell([5, 5, 4])
    print(dir(structure))
    structure.perturb(0.1)

    structure.to(filename="POSCAR_", fmt="poscar")
# 実行
# random()
#
