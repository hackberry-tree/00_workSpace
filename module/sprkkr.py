#!/usr/bin/python
# -*- coding: utf-8 -*-
"""sprkkrのファイルを取り扱う"""
import os
import glob
import copy
import numpy as np
from commopy import Cabinet, Vector

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
            potdat = Pot(pot)
            outdat.data.update(potdat.data)
            if reverse:
                outdat.data.update({'c/a': 1 / potdat.data['c/a']})
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
    """Potファイルを読む (構造の情報)"""
    def __init__(self, fname):
        pot_lines = Cabinet.read_file(fname)
        self.data = self.get_data(pot_lines)

    def get_data(self, pot_lines):
        """
        tetragonalにのみ対応
        c/aとvolumeのみを出力する
        必要あればdictの項目を増やして対応
        その際vaspyのPOSCARを参照したい
        """
        pos_latt = pot_lines.index("LATTICE\n")
        scale = float(pot_lines[pos_latt+4].split()[1])
        cell_lattices = np.array([[float(x) for x in y.split()[1:4]]
                                  for y in pot_lines[pos_latt+5:pos_latt+8]])
        latt_length = [np.linalg.norm(x) * scale for x in cell_lattices]
        cova = latt_length[2] / latt_length[0]
        volume = Vector.get_volume(*cell_lattices) * scale ** 3
        occ = self.get_occupation(pot_lines)
        return({'volume': volume, 'c/a': cova, 'occ': occ})

    def get_occupation(self, pot_lines):
        """
        occupationを読む
        nqを読んで全てのconcをリストするようにしたい(途中)
        現状：一番目のサイトの最後のoccを読む
        """
        nq = pot_lines[9].split() #pylint: disable=C0103
        if nq[0] != 'NQ':
            print('Position of NQ is different')
            exit()
        pos_occ = pot_lines.index("OCCUPATION\n")
        occ = float(pot_lines[pos_occ+2].split()[-1])
        return occ


def main():
    """main"""
    pass

if __name__ == "__main__":
    main()
