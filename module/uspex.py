#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
USPES関連のモジュール
 -Originで構造のツリーを作成
"""
from __future__ import division
import os
import re
import copy
import numpy as np
import random
import pylab
import vaspy
from commopy import Cabinet, Array, Bash, DataBox


def main():
    """main"""
    a = Auxiliary()
    return
    path = '../01_testRun/uspex/'
    path_orig = os.path.join(path, 'origin')
    origin = Origin(path_orig, 1866)
    origin.draw_family_tree(path)


class POSCARS(object):
    """
    UspexのPOSCARSファイルを取り扱う
    元素名がないので、引数にそれを入力
    """
    def __init__(self, path):
        self.path = os.path.dirname(path)
        self.bgp_lines = Cabinet.read_file(path)
        poscar_lines = self.split_structure()
        self.poscars = [{'ID': x['ID'], 'object': vaspy.Poscar(x['lines'])}
                        for x in poscar_lines]

    def expand_poscars(self, elements):
        """
        self.poscarsのPOSCARをディレクトリ(ID名)に展開
        elements名も追記する
        self.path中に作成
        pathを変えたい場合はself.pathを変更
        """
        for poscar in self.poscars:
            Bash.mkdir(os.path.join(self.path, poscar['ID']))
            poscar['object'].elements = elements
            poscar['object'].write_poscar(os.path.join(self.path,
                                                       poscar['ID'], 'POSCAR'))

    def restrict_fractions(self, site_num, fractions):
        """
        ある分率以下の場合、除外する
        Fe-Bの計算に使用
        例として一応残す
        """
        rmvs = []
        for poscar in self.poscars:
            if poscar['object'].get_atom_fractions()[site_num] <= fractions:
                rmvs.append(poscar)
        for rmv in rmvs:
            self.poscars.remove(rmv)

    def remove_elements(self, site):
        """
        POSCAR中のelementが占有するsiteとelement数を消す
        """
        for poscar in self.poscars:
            head_sites = sum(poscar['object'].num_atoms[0: site])
            tail_sites = sum(poscar['object'].num_atoms[0: site+1])
            new_sites = poscar['object'].cell_sites[: head_sites]
            np.append(new_sites, poscar['object'].cell_sites[tail_sites:]) #pylint: disable=E1101
            poscar['object'].cell_sites = new_sites
            poscar['object'].num_atoms.pop(site)





    def split_structure(self):
        """
        各構造毎にPOSCARSを分割
        """
        key = r"^(EA\d+)\s*\n"
        meta = re.compile(key)
        ids = []
        for i in range(0, len(self.bgp_lines)):
            if meta.match(self.bgp_lines[i]):
                num_id = meta.match(self.bgp_lines[i]).group(1)
                ids.append({'ID': num_id, 'pos': i})
        ids.append({'ID': None, 'pos': -1})  # last tail
        poscars = []
        for i in range(0, len(ids) - 1):
            num_id = ids[i]['ID']
            head = ids[i]['pos']
            tail = ids[i+1]['pos']
            lines = self.bgp_lines[head: tail]
            poscars.append({'ID': num_id, 'lines': lines})
        return poscars

class Output(object):
    """
    OUTPUT.txtからparameterを読み取る
    """
    def __init__(self, path):
        self.elements = self.get_elements(path)

    @staticmethod
    def get_elements(path):
        """
        元素名を修得
        """
        lines = Cabinet.read_file(path)
        key = r"\s*There are 3 types of atoms in the system:\s+([\w\s]+)"
        meta = re.compile(key)
        for line in lines:
            if meta.match(line):
                elements = meta.match(line).group(1)
        return elements.split()

class Auxiliary(DataBox):
    """
    USPEXの出力ファイルAuxiliaryFilesを取り扱う
    AuxiliarytyFiles以下のファイル名はデフォルトのまま
    """
    base_dirc = 'AuxiliaryFiles'
    comp_file = 'compositions_nospace.dat'
    energy_file = 'enthalpies_nospace.dat'

    def __init__(self, data):
        DataBox.__init__(self, data)

    @classmethod
    def from_file(cls, path):
        """
        AuxiliaryFile directoryのあるpathを指定
        三つのファイルからデータを読み込む
        """
        return Auxiliary(cls.read_data(path))

    def binary(self):
        """
        binaryのconvex hullをプロット
        """
        self.data.sort(key=lambda x: x['energy'])
        self.data.sort(key=lambda x: x['frac_atoms'][0], reverse=True)
        base0 = self.data[0]['energy']
        self.data.sort(key=lambda x: x['frac_atoms'][1], reverse=True)
        base1 = self.data[0]['energy']
        enthalpy = self['energy'] - (base0 * self['frac_atoms'][:, 0] +
                                     base1 * self['frac_atoms'][:, 1])
        pylab.plot(self['frac_atoms'][:, 1], enthalpy, 'o')
        lines = ""
        i = 0
        for data in self.data:
            lines += " ".join([str(data['frac_atoms'][1]), str(enthalpy[i])]) + "\n"
            i += 1
        print(lines)
        pylab.show()

    def ternary(self):
        """
        ternary
        終端組成の基準energyを0にとったenthalpyを追加する
        """
        self.data.sort(key=lambda x: x['energy'])
        self.data.sort(key=lambda x: x['frac_atoms'][0], reverse=True)
        base0 = self.data[0]['energy']
        self.data.sort(key=lambda x: x['frac_atoms'][1], reverse=True)
        base1 = self.data[0]['energy']
        self.data.sort(key=lambda x: x['frac_atoms'][2], reverse=True)
        base2 = self.data[0]['energy']
        print("base0 is {0} eV/atom".format(base0))
        print("base1 is {0} eV/atom".format(base1))
        print("base2 is {0} eV/atom".format(base2))

        enthalpy = self['energy'] - (base0 * self['frac_atoms'][:, 0] +
                                     base1 * self['frac_atoms'][:, 1] +
                                     base2 * self['frac_atoms'][:, 2])
        self['enthalpy'] = enthalpy

    def separate_bases(self):
        """
        [0, 0, e1], [0, 1, e2], [1, 0, e3]で最小のe*のものを抽出
        そこからのenthalpyに換算してbasesとnot_basesとmeta_stablesをreturnする
        bases: 三角形の頂点3つ
        not_bases: bases以外
        meta_stables: not_basesの中でenthalpyが負のもの
        """
        bases = []
        not_bases = copy.deepcopy(self.data)
        self.data.sort(key=lambda x: x['energy'])
        self.data.sort(key=lambda x: x['frac_atoms'][0], reverse=True)
        bases.append(self.data[0])
        not_bases.remove(self.data[0])
        self.data.sort(key=lambda x: x['frac_atoms'][1], reverse=True)
        bases.append(self.data[0])
        not_bases.remove(self.data[0])
        self.data.sort(key=lambda x: x['frac_atoms'][2], reverse=True)
        bases.append(self.data[0])
        not_bases.remove(self.data[0])
        bases = [[x['frac_atoms'][0],
                  x['frac_atoms'][1],
                  x['enthalpy'] * 96.485344520851] for x in bases]
        not_bases = [[x['frac_atoms'][0],
                      x['frac_atoms'][1],
                      x['enthalpy'] * 96.485344520851] for x in not_bases]
        meta_stables = [x for x in not_bases if x[2] < 0]
        return bases, not_bases, meta_stables

    def __setitem__(self, key, array):
        for i in range(0, len(self.data)):
            if len(self.data) != len(array):
                print("Dimension of data is different")
                break
            try:
                self.data[i][key] = array[i]
            except KeyError:
                self.data[i].update({key: array[i]})

    @classmethod
    def read_data(cls, path):
        """
        AuxiliaryFiles中のデータ収集
        dictのリスト配列として収集する
        """
        comp_path = os.path.join(path, cls.base_dirc, cls.comp_file)
        comp = pylab.loadtxt(comp_path, comments='#')
        total = np.array([[sum(x) for x in comp]])
        frac = comp / total.T
        energy_path = os.path.join(path, cls.base_dirc, cls.energy_file)
        energy = pylab.loadtxt(energy_path, comments="#")
        energy = energy / total[:][0]
        data = [{'frac_atoms': list(x), 'energy': y}
                for x, y in zip(frac, energy)]
        return(data)

    def __getitem__(self, key):
        array = [x[key] for x in self.data]
        return np.array(array)


class Origin(object):
    """
    USPEXの出力ファイルoriginを取り扱う
    python3でしか正しく動作しない
    """
    def __init__(self, fname, descendant):
        """
        attributeは5つ
        individuals: originファイルから読み出した全てのindividual
        descendant: 末裔、これを指定してそれが含まれるlineageを作成
        lineage: descendantから遡ってparentsの情報を収集したもの
        ancestor: lineage中の原種
        has_parents: heredityで進化したindividual
        one_parents: heredity以外のindividual
        """
        self.individuals = self.read_origin_file(fname)
        self.descendant = descendant
        self.lineage = {}
        self.ancestors = []
        self.has_parents = []
        self.one_parents = []
        self.set_lineage()
        average = self.make_family_tree()
        self.make_family_tree(average)

    @staticmethod
    def read_origin_file(fname):
        """
        self.individualsの作成に用いる
        (__init__()で実行)

        Originから読み取ったすべてのindividualの情報を
        dict形式でreturnする
        どの様に進化したか 'background'
        世代 'generation'
        親 'parents'
        子供 'children'(空)を要素に持つ
        またmake_family_tree()で使う'blood'(空)を定義
        """
        individuals = {}
        lines = Cabinet.read_file(fname)
        for line in lines:
            if line.split()[0] == '-' * 7:
                gene = int(line.split()[1].split('generation')[-1])
            else:
                line = line.split(None, 1)
                id_indiv = int(line[0])
                line = line[1].split(',')
                background = line[0].split('\n')[0]
                if background == 'random':
                    parents = None
                else:
                    line = line[1].split()
                    parents = [int(x) for x in line[2:]]
                individuals.update({id_indiv: {'background': background,
                                               'parents': parents,
                                               'generation': gene,
                                               'blood': None,
                                               'children': []}})
        return individuals

    def set_lineage(self):
        """
        末裔 'descendant'を遡る
        set_parents()を繰り返して血族の情報を'self.lineage'にsetする
        'background', 'parents', 'generation'に加えて
        子供 'children'も保管
        """
        descendant = self.descendant
        self.lineage.update({descendant: self.individuals[descendant]})
        self.lineage[descendant].update({'children': None})
        children = [descendant]
        while children:
            self.set_parents(children)
            children, ancestor = self.get_parents_list(children)
            self.ancestors += ancestor
        self.ancestors = sorted(set(self.ancestors), key=self.ancestors.index)

    def set_parents(self, children):
        """
        idから親をたどってself.lineageに親の情報を書き込む
        更に'child'のID情報も加える
        """
        for child in children:
            parents = self.individuals[child]['parents']
            if parents is None:
                pass
            else:
                if len(parents) == 2:
                    self.has_parents.append(child)
                if len(parents) == 1:
                    self.one_parents.append(child)
                for parent in parents:
                    self.lineage.update({parent: {}})
                    self.lineage[parent].update(self.individuals[parent])
                    self.lineage[parent]['children'].append(child)

    def get_parents_list(self, children):
        """
        id_list中のすべての親のリストを作成する
        setを使って重複する親は取り除く
        また親がいないchildrenは
        原種'ansestor'としてreturnする
        """
        parents = []
        ansestor = []
        for child in children:
            add_indiv = self.individuals[child]['parents']
            if add_indiv is None:
                ansestor.append(child)
            else:
                parents += add_indiv
        return sorted(set(parents), key=parents.index), ansestor

    def make_family_tree(self, average=None):
        """
        self.lineage中のindividualに座標を定義して家系図を構成する
        一つの軸は'generation'、もう一軸は血統 'blood'として、
        ancestorsに0~nの整数を与え、子孫の血統はその割合で定義する
        低比率でbloodを混ぜた場合、子孫の中で重なる値が多々出現する
        乱数を使って回避
        点の間隔に疎密が生じるので、bloodの値を再定義する
        (bloodが小さい値から順に等間隔でbloodの値を入れなおす)
        """
        keys = [x for x in self.lineage.keys()]
        blood_list = []
        for i in range(0, len(self.ancestors)):  # ancestorに整数値のbloodを設定
            self.lineage[self.ancestors[i]]['blood'] = i
            blood_list.append(i)
            keys.remove(self.ancestors[i])
        if not average:
            average = 0.5 * (len(self.ancestors) - 1)

        for key in sorted(keys):
            parents = self.lineage[key]['parents']
            if self.lineage[key]['background'] == 'heredity':
                frac = random.random()
                blood_p1 = self.lineage[parents[0]]['blood']
                blood_p2 = self.lineage[parents[1]]['blood']
                blood = (blood_p1 * (0.5 + frac / 10) +
                         blood_p2 * (0.5 - frac / 10))
                blood_list.append(blood)
            elif self.lineage[key]['background'] == 'softmutation':
                blood = self.lineage[parents[0]]['blood']
                if blood < average:
                    blood += 0.01
                else:
                    blood -= 0.01
                blood_list.append(blood)
            else:
                blood = self.lineage[parents[0]]['blood']
            self.lineage[key]['blood'] = blood

        if average == 0.5 * (len(self.ancestors) - 1):
            return self.lineage[self.descendant]['blood']

        blood_list = sorted(set(blood_list))
        blood_redefine = {}

        sum_blood = len(blood_list) - 1
        for i in range(0, len(blood_list)):
            blood_redefine.update({blood_list[i]: i / sum_blood})

        keys = [x for x in self.lineage.keys()]
        for indiv in keys:
            pre_blood = self.lineage[indiv]['blood']
            self.lineage[indiv]['blood'] = blood_redefine[pre_blood]

    def make_family_tree_descendant(self):
        """
        make_family_tree()はancestorsからbloodを決定するが、
        descendantからbloodを決定して行く
        途中が膨らんだtreeになるため不採用にした
        その部分を改善できそうなら使えそう...
        また動作させる為にはlineage中にbranchを定義する必要有り
        しばらく残して、不要そうなら消す
        """
        def set_blood(indiv):
            """
            bloodを
            """
            parents = self.individuals[indiv]['parents']
            if not parents:
                return
            for parent in parents:
                blood = self.lineage[indiv]['blood']
                if self.lineage[parent]['blood'] is None:
                    self.lineage[parent]['blood'] = blood
                branch = self.lineage[indiv]['branch']
                if self.lineage[parent]['branch'] is None:
                    self.lineage[parent]['branch'] = branch

            if len(parents) == 2:
                self.lineage[parents[0]]['blood'] += 1. / branch
                self.lineage[parents[1]]['blood'] -= 1. / branch
                self.lineage[parents[0]]['branch'] += random.random()
                self.lineage[parents[1]]['branch'] += random.random()

        children = [self.descendant]
        while children:
            for child in children:
                set_blood(child)
            children, ancestor = self.get_parents_list(children)
        print(children)
        print(ancestor)

    def draw_family_tree(self, path):
        """
        Draw family tree.
        Make eps file into path
        """
        color_p = {'blue': ('blue', 'cyan'), 'magenta': ('magenta', 'pink'),
                   0: ('#EF6F00', '#FFAF30'), 1: ('#EF7000', '#FFB030'),
                   2: ('#EF4400', '#FF7430'), 3: ('#EF0000', '#FF3030'),
                   4: ('#EF001A', '#FF304A'), 5: ('#EF0044', '#FF3074'),
                   6: ('#EF006F', '#FF309F'), 7: ('#EF009A', '#FF30DA'),
                   8: ('#EF00C4', '#FF30F4'), 9: ('#CF00EF', '#FF30FF'),
                   10: ('#9900EF', '#D930FF'), 11: ('#6F00EF', '#AF30FF'),
                   12: ('#4400EF', '#8430FF'), 13: ('#1A00EF', '#5A30FF'),
                   14: ('#0000EF', '#3030FF'), 15: ('#001AEF', '#305AFF'),
                   16: ('#0045EF', '#3085FF'), 17: ('#006FEF', '#30AFFF'),
                   18: ('#009AEF', '#30DAFF'), 19: ('#00C4EF', '#30F4FF'),
                   20: ('#00EFEF', '#30FFFF'), 21: ('#00EFC4', '#30FFF4'),
                   22: ('#00EF99', '#30FFD9'), 23: ('#00EF6F', '#30FFAF'),
                   24: ('#00EF44', '#30FF84'), 25: ('#00EF1A', '#30FF5A'),
                   26: ('#00EF00', '#30FF30'), 27: ('#45EF00', '#85FF30'),
                   28: ('#6FEF00', '#AFFF30'), 29: ('#9AEF00', '#DAFF30'),
                   30: ('#C4EF00', '#F4FF30'), 31: ('#EFEF00', '#FFFF30'),
                   32: ('#EFC400', '#FFF430'), 33: ('#EF9900', '#FFD930')}

        box = []
        for indiv in self.lineage:
            gene = self.lineage[indiv]['generation']
            blood = self.lineage[indiv]['blood']
            box.append([gene, blood, indiv])
        indiv_points = np.array(box)
        ax1 = pylab.subplot(111)
        ax1.plot(indiv_points[:, 0], indiv_points[:, 1], 'o')

        def plot_heredity(indiv, num_cp):
            """
            Plot heredity with lines
            """
            box = []
            parents = self.lineage[indiv]['parents']
            for key in [parents[0], indiv, parents[1]]:
                gene = self.lineage[key]['generation']
                blood = self.lineage[key]['blood']
                box.append([gene, blood])
            h_lines = np.array(box)
            ax1.plot(h_lines[:, 0], h_lines[:, 1], '-',
                     color=color_p[num_cp][0])
            box = []
            gene = self.lineage[indiv]['generation']
            blood = self.lineage[indiv]['blood']
            box.append([gene, blood])
            h_point = np.array(box)
            ax1.plot(h_point[:, 0], h_point[:, 1], 'o',
                     color=color_p[num_cp][0])

            #ax1.plot(arrayB[:, 0], arrayB[:, 1], 'o')

        i = 0
        for indiv in self.has_parents:
            plot_heredity(indiv, i)
            i += 1
            if i >= 34:
                i = 0

        def connect_indiv(indiv):
            """
            Connect each individuals.
            """
            box = []
            parents = self.lineage[indiv]['parents']
            for key in [parents[0], indiv]:
                gene = self.lineage[key]['generation']
                blood = self.lineage[key]['blood']
                box.append([gene, blood])
            line = np.array(box)
            ax1.plot(line[:, 0], line[:, 1], 'b-')

        for indiv in self.one_parents:
            connect_indiv(indiv)
        pylab.ylim(-.05, 1.05)
        dst_path = os.path.join(path, 'origin.eps')
        pylab.savefig(dst_path)
        pylab.show()


if __name__ == '__main__':
    main()
