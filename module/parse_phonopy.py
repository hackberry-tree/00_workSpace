#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""phonopyのoutputを取り扱う"""
from __future__ import division

import os
import numpy as np

from commopy import DataBox
from vaspy import Poscar

__date__ = "Sep 3 2014"

class Gibbs(DataBox):
    """
    廃止
    gibbs energy vs temperatureのデータを取り扱う
    dict形式{'temp': **, 'gibbs': **}をlistで格納する
    """
    def __init__(self, temperatures, energies):
        if len(temperatures) != len(energies):
            print("Gibbs: data size are different !")
            exit()

        data = []
        for temp, ene in zip(temperatures, energies):
            data.append({'temp': float(temp), 'gibbs': float(ene)})
        DataBox.__init__(self, data)

    def atT(self, temp):
        """
        温度を入力すると線形補間したgibbs energyを算出
        """
        if temp < self.t_min or temp > self.t_max:
            print('Temperature range is out')
            exit()
        i = np.where(self['temp'] == temp)
        if i[0].any():
            return self['gibbs'][i[0][0]]
        else:
            t = self['temp'] - temp
            i = 0
            while t[i] <= 0:
                i += 1
            e0 = self['gibbs'][i-1]
            e1 = self['gibbs'][i]
            t0 = self['temp'][i-1]
            t1 = self['temp'][i]
            de = e1 - e0
            dt = t1 - t0
            gibbs = de / dt * (temp - t0) + e0
            return gibbs

    @property
    def t_min(self):
        """
        温度の最低値を出力
        """
        return min(self['temp'])

    @property
    def t_max(self):
        """
        温度の最高値を出力
        """
        return max(self['temp'])

    @classmethod
    def from_file(cls, fname):
        """
        read from gibbs-temperature.dat file.
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        temperatures = [x.split()[0] for x in lines]
        energies = [x.split()[1] for x in lines]
        return cls(temperatures, energies)

    def alt_gibbs_per_atoms(self, total_atoms):
        """
        単位元素あたりのGibbs energyを追加
        """
        self['gibbs'] /= total_atoms

    # def set_bulk_modulas(self, BulkModulas):
    #     """
    #     BulkModulasをデータに追加
    #     保留
    #     """
    #     self['bulk'] = BulkModulas['bulk']

class BulkModulas(DataBox):
    """
    BulkModulasを取り扱う
    """
    def __init__(self, temperatures, bulk_mods):
        data = []
        for temp, bulk in zip(temperatures, bulk_mods):
            data.append({'temp': float(temp), 'bulk': float(bulk)})
        DataBox.__init__(self, data)

    @classmethod
    def from_file(cls, fname):
        """
        read from bulk_modulus-temperature.dat
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        temperatures = [x.split()[0] for x in lines]
        bulk_mods = [x.split()[1] for x in lines]
        return cls(temperatures, bulk_mods)



class Composition(object):
    """
    組成を取り扱う
    arguments:
        elements: 元素のlist
        num_atoms: 元素の数のlist
    """
    def __init__(self, elements, num_atoms, fract=False):
        self.elements = elements
        self.num_atoms = num_atoms
        self.total_atoms = sum(num_atoms)
        self.fract = fract

    def __len__(self):
        return len(self.dict)

    @property
    def dict(self):
        """
        dict形式で出力
        """
        if self.fract:
            d = self.total_atoms
        else:
            d = 1
        return {x: y / d for x, y in zip(self.elements, self.num_atoms)}

    @classmethod
    def from_poscar(cls, fname):
        """
        POSCARからobjectを作成
        """
        poscar = Poscar(fname)
        return cls(poscar.elements, poscar.num_atoms)

    def __getitem__(self, element):
        """
        元素が含まれている場合はその数をreturn
        無ければ0をreturn
        """
        try:
            return self.dict[element]
        except KeyError:
            return 0

class MurnaghanEoS(DataBox):
    """
    Murnaghan equation of stateのfittingパラメータを取り扱う
    ここからgibbsのfree energyも読めるのでgibbsと統合
    """
    def __init__(self, temps, fit_params):
        data = []
        for temp, param in zip(temps, fit_params):
            [G, B, B1, V0] = param
            data.append({'temp': float(temp), 'G_P0': float(G),
                         'gibbs': float(G),
                         'B': float(B), 'B1': float(B1), 'V0': float(V0)})
        DataBox.__init__(self, data)
        self.is_per_atom = False
        self.pressure = 0

    @classmethod
    def from_file(cls, fname):
        """
        helmholtz-volume.datから読み込む
        """
        with open(fname, 'r') as rfile:
            line = rfile.read()
        temps = []
        fit_params = []
        for item in line.split('Temperature:')[1:]:
            temps.append(float(item.split('\n')[0]))
            fit_params.append(item.split('\n')[1].split()[2:])
        return cls(temps, fit_params)

    def alt_per_atom(self, num_atoms):
        """
        単位元素あたりの値に換算
        """
        num = {True: 1, False: num_atoms}
        self['gibbs'] /= num[self.is_per_atom]
        self['G_P0'] /= num[self.is_per_atom]
        self['V0'] /= num[self.is_per_atom]
        self.is_per_atom = True

    def alt_atP_without_B(self, pressure):
        """
        圧力(GPa)におけるenergyを追加
        """
        self['gibbs'] = self['G_P0'] + pressure / 160.2 * self['V0']

    def alt_atP(self, pressure):
        """
        圧力(GPa)におけるenergyを追加
        """
        P = pressure / 160.2
        B = self['B']
        B1 = self['B1']
        V0 = self['V0']
        V_P = V0 / ((B1 * P / B + 1) ** (1 / B1))
        V0ovV = 1 / ((B1 * P / B + 1) ** (1 / B1))

        dE = B * V_P / (B1 * (B1 - 1)) * (B1 * (1 - V0ovV) + V0ovV ** B1 - 1)

        self['gibbs'] = self['G_P0'] + dE + P * V_P

    @property
    def t_min(self):
        """
        温度の最低値を出力
        """
        return min(self['temp'])

    @property
    def t_max(self):
        """
        温度の最高値を出力
        """
        return max(self['temp'])


    def GatT(self, temp):
        """
        温度を入力すると線形補間したgibbs energyを算出
        """
        if temp < self.t_min or temp > self.t_max:
            print('Temperature range is out')
            exit()
        i = np.where(self['temp'] == temp)
        if i[0].any():
            return self['gibbs'][i[0][0]]
        else:
            t = self['temp'] - temp
            i = 0
            while t[i] <= 0:
                i += 1
            e0 = self['gibbs'][i-1]
            e1 = self['gibbs'][i]
            t0 = self['temp'][i-1]
            t1 = self['temp'][i]
            de = e1 - e0
            dt = t1 - t0
            gibbs = de / dt * (temp - t0) + e0
            return gibbs

class MurnaghanWithComp(object):
    """
    MurnaghanEoSとCompositionを併せたobjectを生成する
    arguments:
        murnaghan: MurnaghanEoS object
        comp: Composition object
        is_normalized: 規格化してあるデータかどうか (bool)
    """
    def __init__(self, murnaghan, comp, is_normalized=False):
        comp.fract = True
        num = {True: 1, False: comp.total_atoms}
        murnaghan.alt_per_atom(num[is_normalized])
        self.murnaghan = murnaghan
        self.comp = comp

    @classmethod
    def from_directory(cls, dirc, poscar='POSCAR',
                      murnaghan_dat='helmholtz-volume.dat'):
        """
        POSCARとhelmholtz-volume.datのあるdirctoryを指定して
        objectを生成
        """
        poscar = os.path.join(dirc, poscar)
        murnaghan_dat = os.path.join(dirc, murnaghan_dat)
        comp = Composition.from_poscar(poscar)
        murnaghan = MurnaghanEoS.from_file(murnaghan_dat)
        return cls(murnaghan, comp)

    def GatT_with_comp(self, temp, elem_list):
        """
        elem_list = [comp_x, comp_y, comp_z]として、
        [comp_x, comp_y, energy]のフォーマットで出力
        convex_hullモジュールに持ち込む
        """
        return [self.comp[elem_list[0]],
                self.comp[elem_list[1]],
                self.murnaghan.GatT(temp)]

class GibbsWithComp(object):
    """
    GibbsとCompositionを併せたobjectを生成する
    arguments:
        gibbs: Gibbs object
        comp: Composition object
        is_normalized: 規格化してあるデータかどうか (bool)
    破棄: murnaghanで代用
    """
    def __init__(self, gibbs, comp, is_normalized=False):
        comp.fract = True
        if is_normalized:
            gibbs.alt_gibbs_per_atoms(1)
            self.gibbs = gibbs
            self.comp = comp
        else:
            gibbs.alt_gibbs_per_atoms(comp.total_atoms)
            self.gibbs = gibbs
            self.comp = comp

    @classmethod
    def from_dirc(cls, dirc, poscar='POSCAR',
                  gibbs_dat='gibbs-temperature.dat'):
        """
        POSCARとgibbs-temperature.datのあるdirctoryを指定して
        objectを生成
        """
        poscar = os.path.join(dirc, poscar)
        gibbs_dat = os.path.join(dirc, gibbs_dat)
        comp = Composition.from_poscar(poscar)
        gibbs = Gibbs.from_file(gibbs_dat)
        return cls(gibbs, comp)

    def atT_with_comp(self, temp, elem_list):
        """
        elem_list = [comp_x, comp_y, comp_z]として、
        [comp_x, comp_y, energy]のフォーマットで出力
        convex_hullモジュールに持ち込む
        """
        return [self.comp[elem_list[0]],
                self.comp[elem_list[1]],
                self.gibbs.atT(temp)]


class References(object):
    """
    sigle elementをkeyにgibss　energyを呼び出せるようにしたobject
    reference energyとして利用
    referenceとしては単体組成のみに対応
    arguments:
        refs: {element: MurnaghanEoS_obj}
    """
    def __init__(self, refs):
        self.elements = list(refs.keys())

        # 温度データが同じであることを確認
        # ToDo: 異なる場合はat_T methodで揃える
        for i in range(1, len(self.elements)):
            if (refs[self.elements[0]]['temp'] != \
                refs[self.elements[i]]['temp']).all():
                print("Temperatures are different !")
                exit()
        self.temp = refs[self.elements[0]]['temp']
        self.refs = refs

    @classmethod
    def from_directories(cls, dir_list):
        """
        dircのlistを指定、dircからobjを生成
        """
        refs = {}
        for dirc in dir_list:
            mwc = MurnaghanWithComp.from_directory(dirc)
            if len(mwc.comp) != 1:
                print('ref is not pure element')
                exit()
            else:
                refs.update({mwc.comp.elements[0]: mwc.murnaghan})
        return cls(refs)

    def at_comp(self, comp, pressure=0):
        """
        referenceを元に組成位置の基準gibbs energyをnp.array形式で出力
        comp: {element: num}のdict形式
        """
        sum_atoms = sum(list(comp.values()))
        ref_at_comp = np.zeros(self.temp.shape)
        for elem in comp:
            self.refs[elem].alt_atP(pressure)
            ref_at_comp += self.refs[elem]['gibbs'] * comp[elem] / sum_atoms
        return ref_at_comp


class Enthalpy(DataBox):
    """
    enthalpyに換算したdata
    arguments:
        compound: 換算元のデータ(MurnaghanWithComp_obj)
        refs: 基準のデータ(References_obj)
    """
    def __init__(self, compound, refs):
        if (compound.murnaghan['temp'] != refs.temp).all():
            print('Enthalpy: Temperatures are different !')
            exit()
        temp = [{'temp': x} for x in refs.temp]
        self.compound = compound
        self.murnaghan = compound.murnaghan
        self.refs = refs
        self.comp = compound.comp
        DataBox.__init__(self, temp)
        self['enthalpy'] = self.enthalpy()

    def enthalpy(self, pressure=0):
        """
        enthalpyをreturn
        """
        self.compound.murnaghan.alt_atP(pressure)
        enthalpy = (self.murnaghan['gibbs'] -
                    self.refs.at_comp(self.comp.dict, pressure))
        return enthalpy

    def atT(self, temp):
        """
        温度を入力すると線形補間したenthalpyを算出
        """
        # if temp < self.t_min or temp > self.t_max:
        #     print('Temperature range is out')
        #     exit()
        i = np.where(self['temp'] == temp)
        if i[0].any():
            return self['enthalpy'][i[0][0]]
        else:
            t = self['temp'] - temp
            i = 0
            while t[i] <= 0:
                i += 1
            e0 = self['enthalpy'][i-1]
            e1 = self['enthalpy'][i]
            t0 = self['temp'][i-1]
            t1 = self['temp'][i]
            de = e1 - e0
            dt = t1 - t0
            enthalpy = de / dt * (temp - t0) + e0
            return enthalpy

    def atT_with_comp(self, temp, elem_list):
        """
        elem_list = [comp_x, comp_y, comp_z]として、
        [comp_x, comp_y, energy]のフォーマットで出力
        convex_hullモジュールに持ち込む
        """
        return [self.comp[elem_list[0]],
                self.comp[elem_list[1]],
                self.atT(temp)]
