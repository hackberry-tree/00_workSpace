#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
input filesを取り扱う
Incar, Posca, Potcar, Kpoints
"""
from __future__ import division
import os
import pymatgen.io.vaspio as Incar
from pymatgen.io.vaspio_set import MITVaspInputSet
from pymatgen.core.structure import Structure

class ProduceInputs(object):
    """
    Inputsファイルの生成を取り纏める
    """
    SAXIS_LIST = [[0, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]]

    def __init__(self, structure, inputset=MITVaspInputSet):
        self.inputset = inputset
        self.structure = structure

    @classmethod
    def from_file(cls, src, inputset=MITVaspInputSet):
        """
        初期構造をfileから読む
        """
        structure = Structure.from_file(src)
        return cls(structure, inputset)

    def write_relax(self, mode, dst='.'):
        """
        構造緩和計算
        mode:
            'all': isif = 3
            'volume': isif = 7
            'ion': isif = 2
        """
        isif = {'all': 3,
                'volume': 7,
                'ion': 2}[mode]
        inputset = self.inputset()
        inputset.kpoints_settings['grid_density'] = 2000
        inputset.incar_settings.pop('LORBIT')
        inputset.incar_settings['NSW'] = 10
        inputset.incar_settings['NELM'] = 60
        inputset.incar_settings['ISIF'] = isif
        inputset.incar_settings['SYSTEM'] = "relax_" + mode
        inputset.write_input(self.structure, dst)

    def write_static(self, dst='.'):
        """
        staticの計算
        """
        inputset = self.inputset()
        inputset.kpoints_settings['grid_density'] = 4000
        inputset.incar_settings.pop('IBRION')
        inputset.incar_settings.pop('ISIF')
        inputset.incar_settings.pop('LORBIT')
        inputset.incar_settings.pop('NSW')
        inputset.incar_settings['ISTART'] = 0
        inputset.incar_settings['ICHARG'] = 2
        inputset.incar_settings['LCHARG'] = False
        inputset.incar_settings['NELM'] = 60
        inputset.incar_settings['NPAR'] = 1
        inputset.incar_settings['SYSTEM'] = "static"
        inputset.write_input(self.structure, dst)

    def write_non_collinear(self, dst='.'):
        """
        mae計算に使用するnon-collinear計算
        pymatgenでnon-collinearのMAGMOM表記が不可なので
        一旦出力後再修正する
        """
        inputset = self.inputset()
        inputset.kpoints_settings['grid_density'] = 4000
        inputset.incar_settings.pop('IBRION')
        inputset.incar_settings.pop('ISIF')
        inputset.incar_settings.pop('LORBIT')
        inputset.incar_settings.pop('NSW')
        inputset.incar_settings.pop('SIGMA')
        inputset.incar_settings['ALGO'] = 'Normal'
        inputset.incar_settings['EDIFF'] = 1e-6
        inputset.incar_settings['GGA_COMPAT'] = False
        inputset.incar_settings['ICHARG'] = 2
        inputset.incar_settings['ISTART'] = 0
        inputset.incar_settings['ISMEAR'] = -5
        inputset.incar_settings['ISYM'] = 0
        inputset.incar_settings['LCHARG'] = True
        inputset.incar_settings['LMAXMIX'] = self._get_lmaxmix()
        inputset.incar_settings['LWAVE'] = True
        inputset.incar_settings['NELM'] = 60
        inputset.incar_settings['NPAR'] = 1



        inputset.incar_settings['SYSTEM'] = "non-collinear"
        inputset.write_input(self.structure, dst)

        incar = inputset.get_incar(self.structure)
        incar['MAGMOM'] = self._magmom_nc_zaxis(incar['MAGMOM'])
        incar.write_file(os.path.join(dst, 'INCAR'))


    @staticmethod
    def _magmom_nc_zaxis(magmom):
        """
        polarized magmomを001方向に向けたnon_collinear表記のmagmomに変換
        """
        return [m for v in zip([0]*len(magmom), [0]*len(magmom), magmom)
                for m in v]

    def _get_lmaxmix(self):
        """
        原子番号からlmaxmixを修得
        """
        if any([el.Z > 56 for el in self.structure.composition]):
            return 6
        elif any([el.Z > 20 for el in self.structure.composition]):
            return 4
        return 2





# class Incar(pmgvaspio.Incar):
#     """
#     Incarを取り扱う
#     """
#     def __init__(self, params=None):
#         """
#         Creates an Incar object.

#         Args:
#             params (dict): A set of input parameters as a dictionary.
#         """
#         super(Incar, self).__init__(params)


#     @classmethod
#     def soc_from_nc(cls, src, dst):
#         """
#         custodianによって修正されたINCARを元にINCARをremakeする
#         INCAR_presoc_ncからINCAR_soc***への変換
#         """
#         diff_soc_nc = ['EDIFF', 'ICHARG', 'ISTART', 'LCHARG', 'LNONCOLLINEAR',
#                        'LSORBIT', 'LWAVE', 'MAGMOM', 'SAXIS']
#         src_inc = cls.from_file(src)
#         dst_inc = cls.from_file(dst)

#         # INCARの差分を取る
#         diff = src_inc.diff(dst_inc)['Different']
#         # socとnc計算由来の差分を除去
#         diff.pop(*diff_soc_nc)
#         diff = {key: values['INCAR2'] for key, values in diff.items()}
#         dst_inc.update(diff)
#         dst_inc.write_file(dst)
