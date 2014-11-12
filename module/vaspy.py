#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module handles VASP Files
inputを作成するIncar, Posca, Potcar, Kpoints
outputを読むOszicar, Outcar
ToDo:INCARのdefault objectのようなものを作成して整理したい
"""
from __future__ import print_function
# from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import generators
# import sys
# if sys.version < '3':
#     text_type = unicode
#     binary_type = str
# else:
#     text_type = str
#     binary_type = bytes

import os
import re
import math
import copy
import socket
import shutil
import numpy as np
import solid
from commopy import Cabinet, Vector, Bash

#================================Gloval Values================================#
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
#================================Gloval Values================================#


def main():
    """
    Execute MakeInputs.all() at current directory.
    """
    MakeInputs.all('.')


class Poscar(object):
    """
    This class manages VASP POSCAR file.
    elements: elements_list
    num_atoms: num of elements_list
    cell_lattices: lattice parameters 3 * 3 array
    cell_deg
    Correspondece only vasp 5
    Noneを引数にいれるとtemplateファイルを読んでobjを作成する
    引数poscarはlist型のlinesにも対応
    """
    def __init__(self, poscar='POSCAR'):
        if type(poscar) is str:
            try:
                poscar_lines = Cabinet.read_file(poscar)
            except IOError:
                print("error: vaspy.Poscar could not "
                      "find '{0}' file !!!".format(poscar))
                exit()
        elif type(poscar) is list:
            poscar_lines = poscar
        elif poscar is None:
            print("POSCAR was not read !!! (Template POSCAR is loaded !!!)")
            poscar = os.path.join(MODULE_DIR,
                                  '../sorce/originalsVASP', 'poscar')
            poscar_lines = Cabinet.read_file(poscar)
        self.poscar_title = poscar_lines[0]
        self.cell_scale = float(poscar_lines[1])
        self.cell_lattices = Cabinet.conv_lines2array(poscar_lines[2:5])
        # self.cell_latticesはarrayとして読み込む
        if poscar_lines[5].split()[0].isdigit():  # vasp4
            self.elements = None
            self.num_atoms = [int(x) for x in poscar_lines[5].split()]
            i = sum(self.num_atoms)
            sites = [[float(x) for x in y.split()[0:3]]
                     for y in poscar_lines[7:7+i]]
            self.cell_sites = np.array(sites)
            self.vasp_version = 4
        else:
            self.elements = poscar_lines[5].split()  # vasp5
            self.num_atoms = [int(x) for x in poscar_lines[6].split()]
            i = sum(self.num_atoms)
            sites = [[float(x) for x in y.split()[0:3]]
                     for y in poscar_lines[8:8+i]]
            self.cell_sites = np.array(sites)
            self.vasp_version = 5

    def __str__(self):
        lines = self.poscar_title
        lines += "  {0:.16f}\n".format(self.cell_scale)
        for latt in self.cell_lattices:
            lines += "  {0[0]:.16f}  {0[1]:.16f}  {0[2]:.16f}\n".format(latt)
        lines += "  {0}\n".format("  ".join(self.elements))
        num_atoms = [str(x) for x in self.num_atoms]
        lines += "  {0}\n".format("  ".join(num_atoms))
        lines += "Direct\n"
        for site in self.cell_sites:
            lines += "  {0[0]:.16f}  {0[1]:.16f}  {0[2]:.16f}\n".format(site)
        return lines

    def get_atom_fractions(self):
        """atomの分率をreturn"""
        sum_atoms = float(sum(self.num_atoms))
        fractions = [x / sum_atoms for x in self.num_atoms]
        return fractions

    def get_cell_volume(self):
        """
        Return cell volume
        """
        volume = Vector.get_volume(*self.cell_lattices) * self.cell_scale ** 3
        return volume

    def get_lattice_length(self):
        """
        Return lattice parameters read from POSCAR
        lattices = [a, b, c]
        """
        lattices = []
        for latt in self.cell_lattices:
            length = np.linalg.norm(latt) * self.cell_scale
            lattices.append(length)
        return lattices

    def get_cell_angle(self):
        """
        return lattices angle of unit cell in degree.
        """
        gamma = (self.cell_lattices[0], self.cell_lattices[1])
        alpha = (self.cell_lattices[1], self.cell_lattices[2])
        beta = (self.cell_lattices[2], self.cell_lattices[0])
        angles = [Vector.get_angle(x, y) for x, y in (alpha, beta, gamma)]
        return angles

    def alt_c_over_a(self, c_over_a):
        """
        This method change C over A with fixed cell_volume and B-axis.
        normalize_latticeからnormalize
        c/aはこのときlen_c/cell_scaleとなる
        """
        self.normalize_lattice()
        len_c = self.get_lattice_length()[2]
        prev = len_c / self.cell_scale
        self.cell_lattices[2] *= c_over_a / prev
        self.cell_scale /= (c_over_a / prev) ** (1./3)
        print("Previous c/a of {} have changed to {}".format(prev, c_over_a))

    def alt_cell_scale(self, scale):
        """
        Alt cell_scale parameter
        """
        self.cell_scale = scale

    def alt_cell_volume(self, volume):
        """
        Alt cell_volume
        """
        ratio = volume / self.get_cell_volume()
        self.cell_scale *= (ratio) ** (1./3.)

    def normalize_lattice(self):
        """
        a軸の記述を1, 0, 0に規格化する
        暫定的に作ったので値のチェックが別途必要
        siteは分率表記である場合、この変換に作用されない
        """
        scale = self.get_lattice_length()[0]
        length = self.get_lattice_length()
        angle = self.get_cell_angle()

        cos_a, cos_b, cos_g = (math.cos(math.radians(x)) for x in angle)
        sin_g = math.sin(math.radians(angle[2]))

        cos_phi = (cos_a - cos_b * cos_g) / sin_g
        cos_theta = ((sin_g ** 2 - cos_a ** 2 - cos_b ** 2 +
                      2 * cos_a * cos_b * cos_g) ** 0.5) / sin_g

        vec_a = [1, 0, 0]
        vec_b = [cos_g, sin_g, 0]
        vec_c = [cos_b, cos_phi, cos_theta]

        self.cell_lattices[0] = vec_a
        self.cell_lattices[1] = vec_b
        self.cell_lattices[2] = vec_c

        self.cell_scale = scale
        self.cell_lattices[0] *= length[0] / scale
        self.cell_lattices[1] *= length[1] / scale
        self.cell_lattices[2] *= length[2] / scale

    def write_poscar(self, poscar='POSCAR'):
        """
        write_poscar(path)
        Make a 'POSCAR' file at 'path'
        """
        Cabinet.reserve_file(poscar)
        Cabinet.write_file(poscar, str(self))


class Potcar(object):
    """
    This class manages potcar files.
    """
    VASP_POT_DIR = os.environ.get('VASP_POTPAW', '')

    def __init__(self, elements=None):
        """
        Make POTCAR to the 'path' directory.
        pot_lines_list
        Get and set the encut, and rwigs values.
        """
        self.psuedo_pot = self.get_psuedo_pot(elements)
        self.potentials_lines = self.read_potcar()

    def read_potcar(self):
        """
        Several POTCAR files lines are loaded based on self.psuedo_pot list.
        """
        path_list = [os.path.join(self.VASP_POT_DIR, x, 'POTCAR')
                     for x in self.psuedo_pot]
        potentials_lines = [Cabinet.read_file(x) for x in path_list]
        return potentials_lines

    @staticmethod
    def get_psuedo_pot(elements):
        """
        Make psuedo_pot list from POT_DICT.
        """
        psuedo_pot = [solid.POT_DICT[x] for x in elements]
        return psuedo_pot

    def read_rwigs(self):
        """
        Read rwigs from POTCAR lines
        """
        rwigs = []
        for p_lines in self.potentials_lines:
            keywords = r"\s*RWIGS\s*=\s*[\d.]+\s*;\s*RWIGS\s*=\s*([\d.]+)\s*.*"
            meta = re.compile(keywords)
            lines_iter = iter(p_lines)
            line = next(lines_iter)
            while meta.match(line) is None:
                line = next(lines_iter)
            match_line = meta.match(line)
            rwigs.append(match_line.group(1))
        return rwigs

    def read_encut(self):
        """
        Read encut from POTCAR lines
        """
        encut_list = []
        for p_lines in self.potentials_lines:
            keywords = r"\s*ENMAX\s*=\s*([\d.]+)\s*;\s*ENMIN\s*=\s*[\d.]+\s*.*"
            meta = re.compile(keywords)
            lines_iter = iter(p_lines)
            line = next(lines_iter)
            while meta.match(line) is None:
                line = next(lines_iter)
            match_line = meta.match(line)
            encut_list.append(float(match_line.group(1)))
        encut = max(encut_list)
        return encut, encut_list

    def write_potcar(self, path, fname='POTCAR'):
        """
        Make a combined single POTCAR file
        """
        fname = os.path.join(path, fname)
        out_lines = [x for y in self.potentials_lines for x in y]
        Cabinet.write_file(fname, out_lines)

    @staticmethod
    def get_composition(fname='./POTCAR'):
        """
        POTCARから元素を読む
        PAW_PBEから始まる行に元素名が記載されているのでそこを読む
        空行は例外処理でpassする
        """
        lines = Cabinet.read_file(fname)
        elements = []
        for line in lines:
            try:
                if line.split()[0] == 'PAW_PBE':
                    elements.append(line.split()[1].split('_')[0])
            except IndexError:
                pass
        return elements


class Kpoints(object):
    """This class manages KPOINTS file."""
    def __init__(self, cell_lattices, dq, parity='odd'):
        """Set self.kpoints"""
        q_vector = [2*math.pi/x for x in cell_lattices]
        kpoints = [int(round(x / dq)) for x in q_vector]
        kpoints_odd = [x + (1 - x % 2) for x in kpoints]
        kpoints_even = [x + (x % 2) for x in kpoints]
        self.kpoints = {'odd': kpoints_odd, False: kpoints,
                        'even': kpoints_even}.get(parity)

    def alt_odd(self):
        """This attr. changes self.kpoints to odd_number"""
        self.kpoints = [x + (1 - x % 2) for x in self.kpoints]

    def alt_size(self, var):
        """Each kpoints times var """
        self.kpoints = [int(x * var) for x in self.kpoints]

    def write_kpoints(self, fname='KPOINTS', mode='M'):
        """Write KPOINTS file using self.kpoints"""
        if mode == 'M':
            mline = "Monkhorst Pack"
        elif mode == 'G':
            mline = "Gamma"
        kp_lines = ("Automatic mesh\n0\n{0}\n"
                    "  {1[0]}  {1[1]}  {1[2]}\n  0.  0.  0.\n"
                    .format(mline, self.kpoints))
        Cabinet.write_file(fname, kp_lines)

    def make_kpoints_band(self, band):
        pass



class IncarReadWriteMixin(object):
    """Read & Write INCAR file methods"""
    @classmethod
    def from_file(cls, fname):
        """
        Read a Incar file, and make a incar_dict.
        """
        lines = Cabinet.read_file(fname)
        incar_dict = {}
        for line in lines:
            if line[0] not in ('#', '\n'):
                para_list = line.split('#')[0].split('!')[0]
                # ^ remove comment_out ^
                para_list = para_list.split()
                key = para_list[0].lower()
                value_list = para_list[2:]
                incar_dict.update({key: value_list})
        cls.__fix_dict(incar_dict)
        return incar_dict

    @classmethod
    def read_incar(cls, fname):
        """
        ToDo: use from_file()
        """
        return cls.from_file(fname)

    @staticmethod
    def __fix_dict(incar_dict):
        """
        読み込んだincar_dictをintやfloat形式に修正
        """
        for key, value in incar_dict.items():
            if len(value) > 1:
                fixed_val = [Cabinet.conv_str(x) for x in value]
            else:
                fixed_val = Cabinet.conv_str(value[0])
            if fixed_val == '.TRUE.':
                fixed_val = True
            elif fixed_val == '.FALSE.':
                fixed_val = False
            incar_dict.update({key: fixed_val})

    def make_incform_all(self):
        """Make INCAR lines into dict."""
        incar_lines = {}
        for key in self.incar_out_list:
            incar_lines.update({key: self.make_incform(key)})
        return incar_lines

    def make_incform(self, key):
        """Change valuables into INCAR format."""
        try:
            var = self.incar_dict[key]
        except KeyError:
            var = None
        if isinstance(var, list):
            var = [str(x) for x in var]
            line = "{0} = {1}\n".format(key.upper(), "  ".join(var))
            return line
        if isinstance(var, bool):
            true = "{0} = {1}\n".format(key.upper(), ".TRUE.")
            false = "{0} = {1}\n".format(key.upper(), ".FALSE.")
            line = {True: true, False: false}.get(var)
            return line
        if isinstance(var, (int, str)):
            line = "{0} = {1}\n".format(key.upper(), var)
            return line
        if isinstance(var, (float)):
            line = "{0} = {1}\n".format(key.upper(), var)
            return line
        if isinstance(var, type(None)):
            line = "\n"
            return line

    def update(self, extra_dict):
        """
        incar_dictをself.updateで要素追加
        __setitem__中に定義したが、
        keyがincar_out_listになければkeyを追加
        """
        for key in extra_dict:
            self[key] = extra_dict[key]

    def __getitem__(self, key):
        return self.incar_dict[key]

    def __setitem__(self, key, var):
        self.incar_dict[key] = var
        if not key in self.incar_out_list:
            self.incar_out_list.append(key)

    def __str__(self):
        """
        fixed_tagをupdateしてINCARのformatでreturn
        """
        self.update(self.cls_fixed_tag)
        self.update(self.fixed_tag)
        lines_dict = self.make_incform_all()
        lines = ""
        for key in self.incar_out_list:
            lines += lines_dict[key]
        while lines.count('\n\n\n') != 0:
            lines = lines.replace('\n\n\n', '\n\n')
        return lines

    def write_incar(self, fname):
        """Write INCAR file"""
        Cabinet.write_file(fname, str(self))


class IncarSwitchTagsMixin(IncarReadWriteMixin):
    """
    INCARのtagをSwichによって切り替える
    またIncar_objを呼ぶ前にtagを追加しておきたい場合、
    cls_add* methodsを利用する
    初期化するときはcls_initializeを使う
    with exit構文で利用できるように整理したい
    """
    @classmethod
    def cls_add_extratag(cls, extra_tag):
        """
        extra_tagを変更
        最初に読み込むのでENCUTなどは適宜INCARで変更される
        """
        cls.cls_extra_tag.update(extra_tag)

    @classmethod
    def cls_add_fixedtag(cls, fixed_tag):
        """
        fixed_tagを変更
        最後に読み込むので作成されるINCARで全て共通になる
        """
        cls.cls_fixed_tag.update(fixed_tag)

    @classmethod
    def cls_initialize(cls):
        """
        classを初期化する
        """
        cls.cls_fixed_tag = {}
        cls.cls_extra_tag = {}

    def switch_istart_lwave(self, read_sw=False, write_sw=False):
        """
        Set istart, icharg, lwave, lcharg.
        -No read and No write (F, F) use cell relaxation, ibzkp calculation.
        -Read and No write (T, F) use for
         spin orbit and band structure calculations.
        -No read and write (F, T) use for
         pre-spin orbit calculations (spin polarized) and pre-band calculation
        """
        read_dict = {True: {'istart': 1, 'icharg': 11},
                     False: {'istart': 0, 'icharg': 2}}
        write_dict = {True: {'lwave': True, 'lcharg': True},
                      False: {'lwave': False, 'lcharg': False}}
        self.update(read_dict[read_sw])
        self.update(write_dict[write_sw])

    def switch_magnetic(self, mag_sw=True):
        """
        磁性計算と非磁性計算の切り替え
        """
        mag_dict = {True: {'ispin': 2}, False: {'ispin': 1, 'magmom': None}}
        self.update(mag_dict[mag_sw])

    def switch_relax_stracture(self, relax_sw=False, isif=3):
        """
        構造緩和計算を行うかどうかの切り替え
        nsw=10とediffg=-0.005がdefault値
        振動してしまってなかなか収束しない
        0.005は結構厳しいのかもしれない...
        directに緩和させる場合はefiffgを整数(default値)にして
        volume依存性から求める場合はediffgを-0.01に変更した

        ibrionを2に、encutを1.3倍に変更する
        """
        relax_dict = {True: {'encut': self.incar_dict['encut'] * 1.3,
                             'ibrion': 2, 'nsw': 10, 'isif': isif,
                             'ediffg': -0.01},
                      False: {'ibrion': None, 'nsw': None, 'isif': None,
                              'ediffg': None}}
        self.update(relax_dict[relax_sw])

    def switch_mae_calc_condition(self, mae_sw=True, lmaxmix=4,
                                  soc_sw=True, saxis=None):
        """
        MAE計算の為のswich
        mae_sw=Trueの場合isymなどを0に指定
        """
        mae_dict = {True: {'gga_compat': False, 'lmaxmix': lmaxmix,
                           'isym': 0, 'ediff': 1.0e-5},
                    False: {'gga_compat': None, 'lmaxmix': None,
                            'isym': None, 'ediff': None}}
        magmom = self.incar_dict['magmom']
        soc_dict = {True: {'lsorbit': True, 'magmom': None, 'saxis': saxis,
                           'ediff': 1.0e-6},
                    False: {'lsorbit': False, 'magmom': magmom, 'saxis': None}}
        self.update(mae_dict[mae_sw])
        self.update(soc_dict[soc_sw])


class IncarLoadPoscarObj(IncarSwitchTagsMixin):
    """
    Correct parameters from Poscar_object (and Potcar).
    """
    cls_extra_tag = {}
    cls_fixed_tag = {}
    incar_out_list = ['system', 's',
                      'npar', 'prec', 'encut', 's',
                      'ispin', 'magmom', 'lsorbit', 'saxis', 's',
                      'gga_compat', 'lmaxmix', 'isym', 's',
                      'nelm', 'nelmin', 'ediff',
                      'ismear', 'sigma', 's',
                      'ibrion', 'nsw', 'isif', 'ediffg', 's',
                      'istart', 'icharg', 'lwave', 'lcharg',
                      'lorbit', 's',
                      'rwigs']

    def __init__(self, poscar_obj):
        self.elements = poscar_obj.elements
        self.num_atoms = poscar_obj.num_atoms
        potcar = Potcar(self.elements)
        self.incar_dict = {}
        self.update({'rwigs': potcar.read_rwigs()})
        self.update({'encut': potcar.read_encut()[0]})
        self.update({'system': self.get_formula()})
        self.update({'magmom': self.make_magmom()})

        self.update({'npar': 1})
        self.update({'prec': 'Accurate'})
        self.update({'ispin': 2})
        self.update(self.cls_extra_tag)
        self.fixed_tag = {}

    def get_formula(self):
        """
        Get system name as a chemical formular.
        """
        formula = ""
        for element, num in zip(self.elements, self.num_atoms):
            if num == 1:
                num = ''
            formula += "{0}{1}".format(element, num)
        return formula

    def make_magmom(self, mag=3):
        """Make magmom from num_atoms."""
        magmom = [mag for y in self.num_atoms for x in range(0, y)]
        #magmom = ["{0}*{1}".format(x, mag) for x in self.num_atoms]
        #return "  ".join(magmom)
        return magmom

class IncarReadPoscar(IncarLoadPoscarObj):
    """
    Correct parameters from Poscar file (and Potcar).
    """
    def __init__(self, poscar='POSCAR'):
        poscar_obj = Poscar(poscar)
        IncarLoadPoscarObj.__init__(self, poscar_obj)


class MakeInputs(object):
    """Make inputs file of series"""
    @classmethod
    def phonon(cls, path, dim):
        """
        phonopy-qhaの計算に必要なファイルを生成
        magmomの行が長くなる場合があるので、別個指定
        """
        poscar_obj = Poscar(os.path.join(path, 'POSCAR'))
        dim_times = dim[0] * dim[1] * dim[2]
        num_atoms = [x * dim_times for x in poscar_obj.num_atoms]
        magmom = "  ".join(["{0}*{1}".format(x, 3) for x in num_atoms])

        incar_obj = IncarReadPoscar(os.path.join(path, 'POSCAR'))
        incar_obj['magmom'] = magmom
        print('Hello!')
        cls.make_potcar_kpoints_for_phonon(path, dim)
        cls.static_for_phonon(path, incar_obj)

        src_dir = os.path.join(MODULE_DIR, '../source/originalsVASP', 'Calc')
        with open(os.path.join(src_dir, 'script-qha.sh'), 'r') as rfile:
            lines = rfile.readlines()

        dim_str = [str(x) for x in dim]
        dim_times = dim[0] * dim[1] * dim[2]
        lines[2] = "dim=\"{0}\"\n".format(" ".join(dim_str))
        lines[5] = "nunit={0}\n".format(dim_times)
        with open(os.path.join(path, 'script-qha.sh'), 'w') as wfile:
            wfile.write("".join(lines))
        os.chmod(os.path.join(path, 'script-qha.sh'), 0o775)

        hostname = socket.gethostname()
        thomson = "thomson.tagen.tohoku.ac.jp"
        if hostname == thomson:
            shutil.copyfile(os.path.join(src_dir, 'calc_thomson.sh'),
                            os.path.join(path, 'calc.sh'))


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
        methods = ['relax', 'cell', 'ion', 'volume', 'volumeE',
                   'presoc', 'presoc_nc',
                   'ibzkp', 'soc',
                   'dos', 'band', 'static']
        for method in methods:
            getattr(cls, method)(path, incar_obj)
        src_dir = os.path.join(MODULE_DIR, '../source/originalsVASP', 'Calc')
        dst_dir = os.path.join(path, 'Calc')
        Bash.copy_dir(src_dir, dst_dir)

    @staticmethod
    def relax(path, base):
        """volume, ion, cell shapes relaxation"""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=3)
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        incar['nelm'] = 40
        fname = os.path.join(path, 'INCAR_relax')
        incar.write_incar(fname)

    @staticmethod
    def cell(path, base):
        """ion & cell shapes relaxation"""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=4)
        incar['encut'] /= 1.3
        incar['isym'] = 0
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        incar['nelm'] = 40
        fname = os.path.join(path, 'INCAR_cell')
        incar.write_incar(fname)

    @staticmethod
    def ion(path, base):
        """only ion relax"""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=2)
        incar['encut'] /= 1.3
        incar['isym'] = 0
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        incar['nelm'] = 40
        fname = os.path.join(path, 'INCAR_ion')
        incar.write_incar(fname)

    @staticmethod
    def volume(path, base):
        """For volume optimize calculation"""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=7)
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        incar.update({'nelm': 40})
        fname = os.path.join(path, 'INCAR_volume')
        incar.write_incar(fname)

    @staticmethod
    def volumeE(path, base): #pylint: disable=C0103
        """For volume optimize calculation"""
        incar = copy.deepcopy(base)
        incar.switch_relax_stracture(relax_sw=True, isif=7)
        del incar.incar_dict['ediffg']
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        incar.update({'nelm': 40})
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
    def static_for_phonon(path, base):
        """
        For static calculation. (no relaxation)
        with ADDGRID = .TRUE.
        """
        incar = copy.deepcopy(base)
        incar.switch_istart_lwave(read_sw=False, write_sw=False)
        incar.switch_mae_calc_condition(mae_sw=False, lmaxmix=4, soc_sw=False)
        incar['isym'] = 0
        incar['addgrid'] = True
        fname = os.path.join(path, 'INCAR')
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

    @staticmethod
    def make_potcar_kpoints_for_phonon(path, dim):
        """
        POTCAR and KPOINTS files are made from SPOSCAR
        """
        poscar = Poscar(os.path.join(path, 'POSCAR'))
        potcar = Potcar(poscar.elements)
        latt = poscar.get_lattice_length()
        latt_spcell = [x * y for x, y in zip(latt, dim)]
        kpoints = Kpoints(latt_spcell, 0.30, 'even')

        potcar.write_potcar(path)
        kpoints.write_kpoints(os.path.join(path, 'KPOINTS'), 'G')


class Oszicar(object):
    """
    OSZICARからのデータ収集
    """
    def __init__(self, fname='OSZICAR'):
        self.results = self.get_results(fname)

    @staticmethod
    def get_3values(line):
        """
        行内からnswとenergyとmagの3つの値をreturnする
        splitした行の長さからmagmomが3次元かスカラーかを判断
        """
        nsw_num = float(line.split()[0])
        energy = float(line.split()[2])
        if len(line.split()) == 10:
            mag = float(line.split()[9])
            mag = math.fabs(mag)
        elif len(line.split()) == 12:
            mag = [float(x) for x in line.split()[9:12]]
            mag = np.linalg.norm(mag)
        else:
            mag = 0
        return nsw_num, energy, mag

    @classmethod
    def get_results(cls, fname='OSZICAR'):
        """
        iterationの回数、nswの数、energy、magをdict形式でreturn
        緩和毎にlistに追加する
        """

        lines = Cabinet.read_file(fname)
        keywords = r"\s*([\d]+)\s+F=\s*([\d\-\.E\+]+)\s+E0=\s+.*\s+"
        meta = re.compile(keywords)
        keywords2 = r"\s*DAV:\s*([\d]+)\s+.*"
        meta2 = re.compile(keywords2)
        results = []
        for i in range(0, len(lines)):
            if meta.match(lines[i]):
                relax_num, energy, mag = cls.get_3values(lines[i])
                j = 1
                while not meta2.match(lines[i-j]):
                    j += 1
                iter_num = lines[i-j].split()[1]
                results.append({'iter_num': iter_num, 'nsw_num': relax_num,
                                'energy': energy, 'mag': mag})
        if not results:
            last_val = lines[-1].split()
            try:
                if math.fabs(float(last_val[3])) > 1e-5:
                    print("{0} is unfinished with error. ".format(fname))
                    return []
            except ValueError:
                    print("{0} is unfinished with error. ".format(fname))
                    return []
            print("{0} is unfinished but converged. "
                  "(val. of mag is false)".format(fname))
            results.append({'iter_num': int(last_val[1]), 'nsw_num': 1,
                            'energy': float(last_val[2]), 'mag': -100})
        return results


class Outcar(object):
    """
    OUTCARからのデータ修得
    """
    def __init__(self, fname='OUTCAR'):
        self.results = self.get_results(fname)

    def get_results(self, fname='OUTCAR'):
        """
        energy, magの値を修得
        """
        lines = Cabinet.read_file(fname)
        energy = self.get_energy(lines)
        mag = self.get_mag(lines)
        elements = self.get_elements(lines)

        results = {'energy': energy, 'mag': mag, 'elements': elements}
        return results

    @staticmethod
    def get_mag(lines):
        """
        magnetic momentを読む
        """
        key = (r"\s*tot\s+[\d\-\.]+\s+[\d\-\.]+"
               r"\s+[\d\-\.]+\s+([\d\-\.]+)\s*.*")
        meta = re.compile(key)
        pos = lines.index(" magnetization (x)\n")
        i = 0
        while not meta.match(lines[pos+i]):
            i += 1
        mag = float(meta.match(lines[pos+i]).group(1))
        mag = math.fabs(mag)
        return mag

    @staticmethod
    def get_energy(lines):
        """
        energyを読む
        """
        key = r"\s*free\s+energy\s+TOTEN\s+=\s+([\d\-\.]+)\s+eV.*"
        meta = re.compile(key)
        energy = [meta.match(x).group(1) for x in lines if meta.match(x)]
        energy = float(energy[-1])
        return energy

    @staticmethod
    def get_elements(lines):
        """
        元素を読む
        """
        key = r"\s+VRHFIN\s+=\s*(.+):\s+.*"
        meta = re.compile(key)
        elements = [meta.match(x).group(1) for x in lines if meta.match(x)]
        return elements


if __name__ == '__main__':
    main()
