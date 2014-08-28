#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
import vaspy

poscar = vaspy.Poscar('POSCAR')
poscar.normalize_lattice()
poscar.write_poscar('POSCAR')
