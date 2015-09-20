#!/usr/bin/python
# -*- coding: utf-8 -*-
"""sprkkrのファイルを取り扱う"""
from __future__ import division
import os
import glob
import copy
import numpy as np
from commopy import Cabinet, Vector
import solid

class CorrectResults(object):
    """
    一連のoutputを収集する
    pathのlistを引数としてlist中の全ての*.potとOUTPUTを読む
    c/a, volumem, energy, magのself.data (array)を作成する
    """
    def __init__(self, path_list, pot=None, output=None):
        """
        initialize
        """
        fnames = {'pot': '*.pot', 'output': 'OUTPUT'}
        if pot:
            fnames['pot'] = pot
        if output:
            fnames['output'] = output
        self.fnames = fnames
        self.path_list = path_list
        self.data = self.get_data(self.path_list)
        self.output = []

    def get_data(self, path_list, reverse=False):
        """
        reverseをTrueにするとc/aをa/cにして読み込む
        (soc100の計算で使用)
        """
        data = []
        for path in path_list:
            output = os.path.join(path, self.fnames['output'])
            pot = glob.glob(os.path.join(path, self.fnames['pot']))[0]
            outdat = Output(output)
            potdat = Pot.from_file(pot)
            outdat.data.update(potdat.data_pre)
            if reverse:
                outdat.data.update({'c/a': 1 / potdat.data_pre['c/a']})
            data.append(outdat.data)
        return data

    def remove_none(self, key):
        """keyのデータに対し空のものを排除する"""
        tmp_data = copy.deepcopy(self.data)
        for data in self.data:
            if data[key] is None:
                tmp_data.remove(data)
        self.data = tmp_data

    def set_mae(self, soc_dir):
        """
        置換するsoc_dirの名前を指定[001, 100]
        path_listは001側のpathを指定しておく
        """
        path_list2 = [x.replace('001','100') for x in self.path_list]
        data2 = self.get_data(path_list2)
        for data100, data001 in zip(data2, self.data):
            try:
                mae = data100['energy'] - data001['energy']
                mae = mae * 10 ** 6 * 13.6058
            except TypeError:
                mae = None
            data001.update({'mae': mae})


    def __getitem__(self, key):
        array = [x[key] for x in self.data]
        return np.array(array)

    def __setitem__(self, key, array):
        for i in range(0, len(self.data)):
            if len(self.data) != len(array):
                print("Dimension of data is different")
                break
            self.data[i].update({key: array[i]})


    def list_no_converged(self):
        """convergencedになっていない結果をlist upする"""
        no_converged = [x for x in self.data if not x['converged']]
        for item in no_converged:
            print("{0[occ]}/{0[c/a]}".format(item))

    def print_table(self):
        """table形式で出力する"""
        lines = "\t".join(self.output) + "\n"
        for data in self.data:
             lines += "\t".join([str(data[x]) for x in self.output]) + "\n"
        print(lines)


class Output(object):
    """結果を読む"""
    def __init__(self, fname):
         output_lines = Cabinet.read_file(fname)
         self.data = self.get_data(output_lines)

    def get_data(self, output_lines):
        if output_lines[-3].find("SCF - cycle converged !!!!!!!!!") != -1:
            energy = float(output_lines[-3].split()[1])
            converged = True
            mag, orbit = [float(x) for x in output_lines[-4].split()[10:12]]
        else:
            energy = None
            mag = None
            orbit = None
            converged = False
        return {'energy': energy, 'mag': mag, 'orbit': orbit,
                'converged': converged}


class Pot(object):
    """
    Potファイルを読む (構造の情報)
    data_pre は暫定的に残して
    data を作成する
    後で統合する
    """
    def __init__(self, data, lines):
        self.lines = lines
        self.data = data

        self.data_pre = {'c/a': self.get_cova_tetragonal(),
                         'volume': self.data['volume'],
                         'occ': self.data['occ'][0][0]['occ']}

    def __str__(self):
        return "".join(self.lines)

    def write_file(self, fname):
        with open(fname, 'w') as wfile:
            wfile.write(str(self))

    @classmethod
    def from_file(cls, fname):
        """
        pot ファイルから object を生成
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        data = cls.get_data(lines)
        return cls(data, lines)

    def get_rws(self, corr_site=1):
        """
        volume から rws を計算
        radii: default の平均イオン半径
        空間群を考慮していないのでユニットセル中のサイト数が実際と違う場合がある
        その場合は corr_site を使って補正する
        """
        radii = []
        for i in range(0, self.data['NQ']):
            radius = 0
            for elem in self.data['occ'][i]:
                radius += (
                    4/3 * np.pi *
                    solid.ELEMENTS[
                        self.data['types'][elem['IT']]['elem']]['Rwigs'] ** 3 *
                    (elem['occ']))
            radii.append(radius)
        volume = self.data['volume']
        scale = volume / sum(radii) / corr_site

        rws = []
        for i in range(0, self.data['NM']):
            _rws = (scale * radii[i] / 4 * 3 / np.pi) ** (1 / 3)
            rws.append(_rws)
        return rws

    def correct_rws(self, corr_site=1):
        """
        dx, rws, rmt の記述を体積に基づいて修正する
        dx = log(rws/R1)/(JWS-1) で与えられる
        R(1)=1e-6, JRMT=721 は 触らない
        """
        rws = self.get_rws(corr_site)
        pos_mesh = self.lines.index("MESH INFORMATION\n")
        for i in range(self.data['NM']):
            var = self.lines[pos_mesh+i+3].split()
            var[2] = "{0:.10f}".format(np.log(rws[i]/1e-6)/(720))
            var[4] = "{0:.10f}".format(rws[i] * 0.85)
            var[6] = "{0:.10f}".format(rws[i])

            self.lines[pos_mesh+i+3] = (
                "    {0[0]}    {0[1]}    {0[2]}    {0[3]}"
                "    {0[4]}  {0[5]}    {0[6]}\n".format(var))

    @classmethod
    def get_data(cls, pot_lines):
        """
        pot_lines から 必要なデータを parse する
        """
        data = {}
        data.update(cls.get_global(pot_lines))
        data.update(cls.get_latt(pot_lines))
        data.update(cls.get_occupation(pot_lines, data['NQ']))
        data.update(cls.get_types(pot_lines, data['NT']))
        return data

    @staticmethod
    def get_latt(pot_lines):
        """
        格子の情報を parse する
        """
        pos_latt = pot_lines.index("LATTICE\n")
        scale = float(pot_lines[pos_latt+4].split()[1])
        cell_lattices = np.array([[float(x) for x in y.split()[1:4]]
                                  for y in pot_lines[pos_latt+5:pos_latt+8]])
        volume = Vector.get_volume(*cell_lattices) * scale ** 3
        return {'volume': volume, 'latt': cell_lattices, 'latt_scale': scale}

    @staticmethod
    def get_global(pot_lines):
        """
        GLOBAL SYSTEM PARAMETER を parse する
        """
        pos_gl = pot_lines.index("GLOBAL SYSTEM PARAMETER\n")
        nq = int(pot_lines[pos_gl+1].split()[1])
        nt = int(pot_lines[pos_gl+2].split()[1])
        nm = int(pot_lines[pos_gl+3].split()[1])
        return {'NQ': nq, 'NT': nt, 'NM': nm}

    @staticmethod
    def get_occupation(pot_lines, nq):
        """
        OCCUPATION を parse する
        occ = [site1, site2, site3, .. site_nq]
        site* = [{IT: * occ: *}, {IT: *, occ: *}]
        のフォーマットで return する
        """
        pos_occ = pot_lines.index("OCCUPATION\n")
        occ = []
        for i in range(nq):
            occ_var = pot_lines[pos_occ+2+i].split()
            site = [{'IT': int(occ_var[4+j*2]), 'occ': float(occ_var[5+j*2])}
                    for j in range(int(occ_var[3]))]
            occ.append(site)
        return {'occ': occ}

    @staticmethod
    def get_types(pot_lines, nt):
        """
        TYPES を parse する
        ZT の値から {IT: {ZT, elem} を return
        """
        pos_types = pot_lines.index("TYPES\n")
        types = {}
        rev_ele = {v['Z']: k for k, v in solid.ELEMENTS.items()}
        for i in range(nt):
            types_var = pot_lines[pos_types+2+i].split()
            types.update({i+1: {'Z': int(types_var[2]),
                                'elem': rev_ele[int(types_var[2])]}})
        return {'types': types}

    def get_cova_tetragonal(self):
        """
        c/a を return する
        tetragonalにのみ対応
        """

        latt_length = [np.linalg.norm(x) * self.data['latt_scale']
                       for x in self.data['latt']]
        return latt_length[2] / latt_length[0]



def main():
    """main"""
    pass

if __name__ == "__main__":
    main()
