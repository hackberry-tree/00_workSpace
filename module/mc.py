#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Cluster密度を計算
"""
from __future__ import division
import os
import glob
import random
import argparse

import itx
import collect_vasp

class QuadSite(object):
    """
    4点のサイトの情報を記述する
    spins
    sites
    """
    def __init__(self, spins):
        self.spins = CyclicList(spins)
        site0 = [0, 0, 0]
        site1 = [0, 1/2, 1/2]
        site2 = [1/2, 0, 1/2]
        site3 = [1/2, 1/2, 0]
        self.sites = CyclicList([site0, site1, site2, site3])

    @classmethod
    def random(cls, conc):
        spins = [int(2*round(random.uniform(conc/2., 0.5+conc/2.))-1)
                for i in range(4)]
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


class Cell(object):
    """
    周期境界条件を課した3次元のlist
    セルの大きさはsize^3
    QuadSiteのどのsiteを原点にするかをoriginで指定する
    """
    def __init__(self, size, random=True, origin=0):
        self.size = size
        tet = QuadSite.random
        self.cell = CyclicList(CyclicList(CyclicList(tet(0.5)
                                                     for z in range(size))
                                          for y in range(size))
                               for x in range(size))
        self.origin = origin

    def get_tetra_xi(self, coord):
        """
        座標周りの8つのtetragonalクラスターの相関関数を計算する
        """
        s0 = self[coord]
        print(s0)

        # s1 =
        # s2
        # s3

    def from_origin(self, site_id, *pos):
        """
        QuadSiteの原点を番号で指定
        その原点からの相対位置で各spinを参照する
        """
        trans = [self._orig0, self._orig1,
                 self._orig2, self._orig3][site_id](pos[0])
        rpos = [i + j for i, j in zip(pos[0], trans)]
        return self[rpos]

    def _from_origin(self, *pos):
        """
        QuadSiteの原点を番号で指定
        その原点からの相対位置で各spinを参照する
        """
        trans = [self._orig0, self._orig1,
                 self._orig2, self._orig3][self.origin](pos[0])
        rpos = [i + j for i, j in zip(pos[0], trans)]
        return rpos


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
                 2:[0, 0, 1, 1], 3:[0, 1, 0, -1]}[pos[3]%4]
        return trans

    @staticmethod
    def _orig2(pos):
        """
        並進対称ベクトル site2を原点にした場合
        """
        trans = {0: [0, 0, 0, 2], 1: [0, 0, 1, 2],
                 2:[1, 0, 1, -2], 3:[1, 0, 0, -2]}[pos[3]%4]
        return trans

    @staticmethod
    def _orig3(pos):
        """
        並進対称ベクトル site3を原点にした場合
        """
        trans = {0: [0, 0, 0, 3], 1: [0, 1, 0, 1],
                 2:[1, 0, 0, -1], 3:[1, 1, 0, -3]}[pos[3]%4]
        return trans

    def get_conc(self, i):
        """
        A元素の格子中の濃度を算出
        """
        spin = [s for x in self.cell for y in x for z in y for s in z.spins]
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

l = 10**5
print(l)
conc = 0
j = []
for i in range(l):
    j.append(int(2*round(random.uniform(conc/2., 0.5+conc/2.))-1))
print(j[0])
print(j.count(1)/float(l))
rand = QuadSite.random(0.5)
print(rand)

cell = Cell(24)
print(cell.get_conc(1))
print(cell[1,1,1,:] == cell[25,25,25,:])
print(cell.from_origin(0, [1,1,1,1]))
cell.get_tetra_xi([0,0,1])

# test
print(cell.from_origin(1, [0, 0, 0, 0]) == cell.from_origin(2, [-1, 0, 0, 3]))
print(cell.from_origin(2, [0, 0, 0, 3]) == cell.from_origin(3, [0, 0, 0, 2]))
print(cell.from_origin(3, [0, 0, 0, 1]) == cell.from_origin(0, [0, 1, 0, 2]))
print(cell.from_origin(1, [1, 2, 3, 1]) == cell.from_origin(0, [1, 3, 4, 0]))
print(cell.from_origin(2, [0, 0, 0, 0]) == cell.from_origin(0, [0, 0, 0, 2]))
