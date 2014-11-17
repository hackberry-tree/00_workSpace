#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Cluster密度を計算
"""
from __future__ import division
import numpy as np
import math
import random
#import pyximport; pyximport.install(pyimport=True)
from correlationfunction import Tetrahedron

def main():
    cell = CellFCC(10, mode='L10')
    print('concentration')
    print(cell.get_conc(1))
    print('local xi at 0 0 0')
    print(cell.get_tetra_xi([0, 0, 0], 0))
    print('gloval xi')
    print(cell.get_tetra_xi_glob())

class QuadSite(object):
    """
    4点のサイトの情報を記述する
    spins
    sites
    """
    def __init__(self, spins):
        self.spins = spins
        site0 = [0, 0, 0]
        site1 = [0, 1/2, 1/2]
        site2 = [1/2, 0, 1/2]
        site3 = [1/2, 1/2, 0]
        self.sites = [site0, site1, site2, site3]

    @classmethod
    def random(cls, conc):
        spins = [int(2*round(random.uniform(conc/2., 0.5+conc/2.))-1)
                for i in range(4)]
        return cls(spins)

    @classmethod
    def L10(cls, conc=0.5):
        spins = [-1, -1, 1, 1]
        return cls(spins)

    @classmethod
    def L12_A(cls, conc=3/4):
        spins = [-1, 1, 1, 1]
        return cls(spins)

    @classmethod
    def L12_B(cls, conc=1/4):
        spins = [-1, -1, -1, 1]
        return cls(spins)


    def __str__(self):
        return str(self.spins)

    def __getitem__(self, index):
        return self.spins[index]

class CyclicList(list):
    """
    周期境界条件を課したlist
    """
    def __init__(self, item):
        super(CyclicList, self).__init__(item)

    def __getitem__(self, index):
        """
        index%len(self)で周期性を与える
        sliceが使えるようTypeErrorの場合はpassする
        """
        try:
            index = index % len(self)
        except TypeError:
            pass
        return super(CyclicList, self).__getitem__(index)

    def __setitem__(self, index, item):
        """
        index%len(self)で周期性を与える
        """
        try:
            index = index % len(self)
        except TypeError:
            pass
        return super(CyclicList, self).__setitem__(index, item)


class CellFCC(object):
    """
    周期境界条件を課した3次元のlist
    セルの大きさはsize^3
    QuadSiteのどのsiteを原点にするかをoriginで指定する
    """
    # 原点からみた最隣接原子の座標
    NEAREST = [[0, 0, 0, 1], [0, 0, 0, 2], [0, -1, 0, 1], [-1, 0, 0, 2],
               [0, 0, 0, 3], [0, -1, 0, 3], [-1, -1, 0, 3], [-1, 0, 0, 3],
               [0, 0, -1, 1], [0, 0, -1, 2], [0, -1, -1, 1], [-1, 0, -1, 2]]
    # memo 上記の実座標 簡単のため2倍して表記
    COORDINATE_NEAREST = [[0, 1, 1], [1, 0, 1], [0, -1, 1], [-1, 0, 1],
                          [1, 1, 0], [1, -1, 0], [-1, -1, 0], [-1, 1, 0],
                          [0, 1, -1], [1, 0, -1], [0, -1, -1], [-1, 0, -1]]
    # NEARESTを3つ選んで作る8つの四面体
    TETRA_POS = [[0, 1, 4], [3, 0, 7], [2, 3, 6], [1, 2, 5],
                 [8, 9, 4], [11, 8, 7], [10, 11, 6], [9 ,10, 5]]

    # 原点から見た第二隣接原子の座標
    NEXT_N = [[1, 0, 0], [-1, 0, 0],
              [0, 1, 0], [0, -1, 0],
              [0, 0, 1], [0, 0, -1]]
    # 第二隣接一つと4つの最隣接からなる八面体
    # OCTAの順番で第二隣接位置はNEXT_Nの順番と対応
    # 例:OCTA[1]の第二隣接位置はNEXT_N[1](=[-1, 0, 0])
    OCTA_POS = [[1, 9, 4, 5], [3, 11, 6, 7],
                [0, 8, 4, 7], [2, 10, 5, 6],
                [0, 2, 1, 3], [8, 10, 9, 11]]

    def __init__(self, size, mode='random', origin=0):
        self.size = size
        tet = {'random': QuadSite.random, 'L10': QuadSite.L10,
               'L12_A': QuadSite.L12_A, 'L12_B': QuadSite.L12_B,}[mode]
        self.cell = np.array([[[tet(0.5).spins for z in range(size)]
                               for y in range(size)] for x in range(size)])
        self.origin = origin

        self.TETRA = [[self.NEAREST[y] for y in x] for x in self.TETRA_POS]

    def get_tetra_xi(self, coord, site_id=0):
        """
        座標周りの8つのtetragonalクラスターの相関関数を計算する
        coordは(x, y, z)の三次元座標 site_idを別途指定する
        """
        c = coord + [0]  # 末尾に0番を追加
        spin0 = self.from_origin(c, site_id)

        # coord周りのtetrahedronの座標
        tetra = [[[x+y for x, y in zip(c, s)] for s in t] for t in self.TETRA]
        # 四面体クラスターのspin * 8つ
        spins = [[spin0] + self.get_arrange_spins(t, site_id) for t in tetra]
        return Tetrahedron.average(spins)

    def get_tetra_xi_glob(self):
        """
        系内全域の相関関数をreturn
        """
        coords = [[x, y, z] for x in range(-1, self.size-1)
                  for y in range(-1, self.size-1)
                  for z in range(-1, self.size-1)]
        site_ids = [0, 1, 2, 3]
        xis = [self.get_tetra_xi(c, i) for c in coords for i in site_ids]
        return [math.fsum(i)/len(i) for i in zip(*xis)]

    def get_arrange_spins(self, pos_list, site_id=0):
        """
        pos_list中の各座標におけるspinを参照し、spin配列のlistとしてreturnする
        site_idでどの原点から見た座標なのかを指定する
        """
        return [self.from_origin(x, site_id) for x in pos_list]

    def from_origin(self, pos, site_id=0):
        """
        QuadSiteの原点を番号で指定
        その原点からの相対位置で各spinを参照する
        スライスは使えない
        QuadSiteごと組み替える必要があって原理的に難しい
        """
        trans = [self._orig0, self._orig1,
                 self._orig2, self._orig3][site_id](pos)
        rpos = [i + j for i, j in zip(pos, trans)]
        return self[rpos]

    @staticmethod
    def _orig0(_):
        """
        並進対称ベクトル 変換無し
        """
        return [0, 0, 0, 0]

    @staticmethod
    def _orig1(pos):
        """
        並進対称ベクトル site1を原点にした場合
        """
        trans = {0: [0, 0, 0, 1], 1: [0, 1, 1, -1],
                 2: [0, 0, 1, 1], 3: [0, 1, 0, -1]}[pos[3]%4]
        return trans

    @staticmethod
    def _orig2(pos):
        """
        並進対称ベクトル site2を原点にした場合
        """
        trans = {0: [0, 0, 0, 2], 1: [0, 0, 1, 2],
                 2: [1, 0, 1, -2], 3: [1, 0, 0, -2]}[pos[3]%4]
        return trans

    @staticmethod
    def _orig3(pos):
        """
        並進対称ベクトル site3を原点にした場合
        """
        trans = {0: [0, 0, 0, 3], 1: [0, 1, 0, 1],
                 2: [1, 0, 0, -1], 3: [1, 1, 0, -3]}[pos[3]%4]
        return trans

    def get_conc(self, i):
        """
        A元素の格子中の濃度を算出
        """
        spin = [s for x in self.cell for y in x for z in y for s in z]
        return spin.count(i) / len(spin)

    def slice_x(self, x):
        """
        座標xの面でスライスする
        """
        pass

    def slice_y(self, y):
        """
        座標xの面でスライスする
        """
        pass

    def slice_z(self, z):
        """
        座標xの面でスライスする
        """
        pass

    def __getitem__(self, *pos):
        pos = pos[0]
        try:
            return self.cell[pos[0]][pos[1]][pos[2]][pos[3]]
        except IndexError:
            return self.cell[pos[0]][pos[1]][pos[2]]

if __name__ == '__main__':
    main()


