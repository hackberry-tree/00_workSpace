#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Poscarの情報を修得
"""
import vaspy

poscar = vaspy.Poscar('POSCAR')
print("cell volume is {0}".format(poscar.get_cell_volume()))
