#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
core
基礎的なclass集
"""
import math
import numpy as np
from collections import defaultdict

class UnitCell(object):
    """
    unit cell object
    """
    def __init__(self, trans_vect, sites):
        self.trans_vect = trans_vect
        self.sites = sites

    @property
    def volume(self):
        """
        unit cellのvolumeをreturn
        """
        return self.trans_vect.volume

    @property
    def degree(self):
        """
        unit cellの角度をreturn
        """
        return self.trans_vect.angle


class TranslationVector(object):
    """
    unit cellの外枠を記述する
    またvolumeも算出
    augments:
        lattices: [latt_x, latt_y, latt_z, ...] (list of Lattice objects)
    """
    def __init__(self, lattices):
        self.lattices = lattices

    @property
    def dimension(self):
        """dimension of unit cell"""
        return [len(self.lattices), self.lattices[0].dimension]

    @property
    def volume(self):
        """
        体積を算出
        """
        matrix = np.array([x.vector for x in self.lattices])
        return math.fabs(np.linalg.det(matrix))

    @property
    def angle(self):
        """
        3軸の角度を算出
        """
        alpha = self._get_angle(self.lattices[1], self.lattices[2])
        beta = self._get_angle(self.lattices[2], self.lattices[1])
        gamma = self._get_angle(self.lattices[0], self.lattices[1])
        return (alpha, beta, gamma)

    @staticmethod
    def _get_angle(latt_1, latt_2):
        """
        Calculate angle of between two vectors.
        unit: degree
        """
        times_12 = latt_1.length * latt_2.length
        cosine = np.dot(latt_1.vector, latt_2.vector) / times_12
        degree = math.acos(cosine) / (2 * math.pi) * 360
        return degree

    @classmethod
    def from_list(cls, latt_list):
        """
        listからobjectを生成
        """
        x = Lattice(latt_list[0], 'x')
        y = Lattice(latt_list[0], 'y')
        z = Lattice(latt_list[0], 'z')
        return cls([x, y, z])

class Lattice(object):
    """
    格子のベクトルを記述する
    arguments:
        vector: lattice of vector (list of float)
        axis: label for axis (strings)
    """
    def __init__(self, vector, axis=None):
        self.vector = vector
        self.axis = axis

    def __len__(self):
        return len(self.vector)

    @property
    def length(self):
        """lattice length"""
        return np.linalg.norm(np.array(self.vector))

    @property
    def dimension(self):
        """dimension of vector"""
        len(self)

    @property
    def to_dict(self):
        """dict形式で出力"""
        return {'axis': self.axis, 'vector': self.vector}

class Sites(object):
    """
    サイトを記述する
    arguments:
        positions: サイト位置を記したvector
        elements: サイト位置に配置される元素
        positionsとelementsは同じ大きさのlistである必要がある
        elementsを省略した場合、Noneをkeyとするdictを作成
    """
    def __init__(self, positions, elements=None):
        self.sites = defaultdict(list)
        if elements:
            if len(positions) != len(elements):
                print("two is different !")
                exit()
            for pos, elem in zip(positions, elements):
                self.sites[elem] += [pos]
        else:
            self.sites[None] = positions

    def __str__(self):
        lines = ""
        for elem in self.sites:
            for site in self.sites[elem]:
                lines += elem + "\t"
                str_site = [str(x) for x in site]
                lines += "\t".join(str_site) + "\n"
        return lines
