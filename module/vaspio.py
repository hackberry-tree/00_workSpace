#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module handles VASP Files
inputを作成するIncar, Posca, Potcar, Kpoints
outputを読むOszicar, Outcar
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import generators

import os
import re
import math
import copy
import numpy as np
import solid
from commopy import Cabinet, Vector, Bash
from pymatgen.io.vaspio import vasp_input

class Poscar(vasp_input.Poscar):
    def __init__(self, structure, comment=None, selective_dynamics=None,
                 true_names=True, velocities=None, predictor_corrector=None):
        vasp_input.Poscar.__init__(self, structure, comment=None,
                                   selective_dynamics=None, true_names=True,
                                   velocities=None, predictor_corrector=None)

class ProduceVaspInputs(object):
    """Make inputs file of series"""
    @classmethod
    def all(cls, path, incar_obj=None, kp_rx=0.15, kp_soc=0.11):
        """
        All pattern of INCAR files are prepared.
        それぞれのパラメータの変更は
        class method　"cls_add_fixedtag"や"cls_add_extratag"を使うか
        incar_objに書き替えたincarを入力
        省略した場合、path中のPOSCARをbaseに作成
        """
        if not incar_obj:
            incar_obj = IncarReadPoscar(os.path.join(path, 'POSCAR'))
        cls.make_potcar_kpoints(path, kp_rx, kp_soc)
        methods = ['relax', 'cell', 'volume', 'volumeE', 'presoc', 'presoc_nc',
                   'ibzkp', 'soc',
                   'dos', 'band', 'static']
        for method in methods:
            getattr(cls, method)(path, incar_obj)
        src_dir = os.path.join(MODULE_DIR, '../sorce/originalsVASP', 'Calc')
        dst_dir = os.path.join(path, 'Calc')
        Bash.copy_dir(src_dir, dst_dir)

    @staticmethod
    def relax(path, base):
        """For cell optimize calculation"""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=3)
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        fname = os.path.join(path, 'INCAR_relax')
        incar.write_incar(fname)

    @staticmethod
    def cell(path, base):
        """For cell optimize calculation"""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=4)
        incar['encut'] /= 1.3
        incar['isym'] = 0
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        fname = os.path.join(path, 'INCAR_cell')
        incar.write_incar(fname)

    @staticmethod
    def volume(path, base):
        """For volume optimize calculation"""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=7)
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        fname = os.path.join(path, 'INCAR_volume')
        incar.write_incar(fname)

    @staticmethod
    def volumeE(path, base): #pylint: disable=C0103
        """For volume optimize calculation"""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=7)
        del incar.incar_dict['ediffg']
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        fname = os.path.join(path, 'INCAR_volumeE')
        incar.write_incar(fname)

    @staticmethod
    def cell_nonmag(path, base):
        """For cell optimize calculation with nonmag ver."""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=3)
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        incar.switch_magnetic(False)
        fname = os.path.join(path, 'INCAR_cell_nonmag')
        incar.write_incar(fname)

    @staticmethod
    def volume_nonmag(path, base):
        """For volume optimize calculation with nonmag ver."""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=7)
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        incar.switch_magnetic(False)
        fname = os.path.join(path, 'INCAR_volume_nonmag')
        incar.write_incar(fname)

    @staticmethod
    def presoc(path, base):
        """
        For spin-polarized calculation to generate WAVECAR and CHGCAR
        for non-self consistent soc calculations.
        """
        incar = copy.deepcopy(base)
        incar.switch_istart_lwave(read_sw=False, write_sw=True)
        incar.switch_mae_calc_condition(mae_sw=True, lmaxmix=4, soc_sw=False)
        incar.incar_dict.update({'nelm': 150})
        incar.incar_dict.update({'ismear': -5})
        fname = os.path.join(path, 'INCAR_presoc')
        incar.write_incar(fname)

    @staticmethod
    def ibzkp(path, base):
        """
        Make IBZKP file from 1 iteration calculation of soc condition.
        """
        incar = copy.deepcopy(base)
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        incar.switch_mae_calc_condition(mae_sw=True, lmaxmix=4,
                                        soc_sw=True, saxis=[1, 0, 0])
        incar.incar_dict.update({'nelm': 1})
        incar.incar_dict.update({'nelmin': 1})
        incar.incar_dict.update({'ismear': -5})
        fname = os.path.join(path, 'INCAR_ibzkp')
        incar.write_incar(fname)

    @staticmethod
    def presoc_nc(path, base):
        """
        Noncollinear caluculation to generate WACECAR and CHGCAR
        for non-self consistent soc calculation.
        見やすくするためmagmomはstrで作成
        """
        incar = copy.deepcopy(base)
        incar.switch_istart_lwave(read_sw=False, write_sw=True)
        incar.switch_mae_calc_condition(mae_sw=True, lmaxmix=4, soc_sw=False)
        incar.incar_dict.update({'lnoncollinear': True})
        magmom_tmp = [['0', '0', str(x)] for x in incar.incar_dict['magmom']]
        magmom_3d = ""
        for mag in magmom_tmp:
            magmom_3d += " ".join(mag) + "  "
        incar.incar_dict.update({'magmom': magmom_3d})
        incar.incar_dict.update({'nelm': 150})
        incar.incar_dict.update({'ismear': -5})
        incar.incar_out_list.append('lnoncollinear')
        fname = os.path.join(path, 'INCAR_presoc_nc')
        incar.write_incar(fname)

    @staticmethod
    def soc(path, base):
        """
        For soc calculation.
        """
        incar = copy.deepcopy(base)
        incar.switch_istart_lwave(read_sw=True, write_sw=False)
        incar.incar_dict.update({'ismear': -5})
        saxis_list = ['001', '100', '110', '111']
        for saxis in saxis_list:
            direction = [int(x) for x in saxis]
            incar.switch_mae_calc_condition(mae_sw=True, lmaxmix=4,
                                            soc_sw=True, saxis=direction)
            fname = os.path.join(path, 'INCAR_soc{0}'.format(saxis))
            incar.write_incar(fname)

    @staticmethod
    def dos(path, base):
        """
        For pre-band calculation.
        ISMEAR shoud be 1. (maybe)
        """
        incar = copy.deepcopy(base)
        incar.switch_istart_lwave(read_sw=False, write_sw=True)
        incar.switch_mae_calc_condition(mae_sw=True, lmaxmix=4, soc_sw=False)
        incar.incar_dict.update({'nelm': 150})
        incar.incar_dict.update({'lorbit': 2})
        incar.incar_dict.update({'ismear': 1})
        incar.incar_dict.update({'sigma': 0.02})
        fname = os.path.join(path, 'INCAR_dos')
        incar.write_incar(fname)

    @staticmethod
    def band(path, base):
        """
        For band calculation.
        ISMEAR shoud be 1. (maybe)
        """
        incar = copy.deepcopy(base)
        incar.switch_istart_lwave(read_sw=True, write_sw=False)
        incar.switch_mae_calc_condition(mae_sw=True, lmaxmix=4, soc_sw=False)
        incar.incar_dict.update({'nelm': 150})
        incar.incar_dict.update({'lorbit': 2})
        incar.incar_dict.update({'ismear': 1})
        incar.incar_dict.update({'sigma': 0.02})
        incar.incar_dict.update({'ediff': 1.0e-6})
        fname = os.path.join(path, 'INCAR_band')
        incar.write_incar(fname)

    @staticmethod
    def static(path, base):
        """
        For static calculation. (no relaxation)
        """
        incar = copy.deepcopy(base)
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        incar.switch_mae_calc_condition(mae_sw=False, lmaxmix=4, soc_sw=False)
        incar['isym'] = 0
        fname = os.path.join(path, 'INCAR_static')
        incar.incar_dict.update({'encut': 400})
        incar.write_incar(fname)

    @staticmethod
    def make_potcar_kpoints(path, relax=0.15, soc=0.11):
        """
        POTCAR and KPOINTS files are made from POSCAR condition.
        """
        poscar = Poscar(os.path.join(path, 'POSCAR'))
        potcar = Potcar(poscar.elements)
        kpoints_relax = Kpoints(poscar.get_lattice_length(), relax)
        kpoints_reduc = Kpoints(poscar.get_lattice_length(), relax*2)
        kpoints_soc = Kpoints(poscar.get_lattice_length(), soc)
        potcar.write_potcar(path)
        kpoints_relax.write_kpoints(os.path.join(path, 'KPOINTS_relax'))
        kpoints_reduc.write_kpoints(os.path.join(path, 'KPOINTS_relax_reduced'))
        kpoints_soc.write_kpoints(os.path.join(path, 'KPOINTS_soc'))
