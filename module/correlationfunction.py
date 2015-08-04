#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
相関関数を計算する
要素はnumpyのarrayで
"""
import math
import numpy as np
from itertools import combinations, chain
from makeSeries import MakePattern

class SpinOperator(object):
    """
    spinの計算のmethod
    """
    @classmethod
    def sum_spins(cls, spins):
        """
        spinの積をとってsumをreturnする
        """
        return np.array([cls.times_spins(s) for s in spins]).mean()

    @staticmethod
    def times_spins(spins):
        """
        全てのspinの積をretrun
        """
        # i = 1
        # for spin in spins:
        #     i *= spin
        # return i
        spins = np.array(spins)
        return spins.prod()


class Octahedron(SpinOperator):
    """
    四面体と異なってサイトの位置関係で相関関数が異なる
    spinsは(0,1) (2,3) (4,5)が第二隣接の組み合わせとして入力する
    spinsを第二隣接同士の各組み合わせで平均をとって[-1/0/1]*3のlistに
    変換して処理した方が早いかもしれない
    パターンとして3^3=27通りを分類分けすることになる
    """
    def __init__(self, ecis):
        self.ecis = ecis

    @staticmethod
    def _nn_combi(spins):
        """
        第二隣接同士で組み合わせたlistをreturnする
        """
        return [spins[0:2],
                spins[2:4],
                spins[4:6]]

    @staticmethod
    def _null(_):
        """
        nullクラスター
        """
        return 1

    @classmethod
    def _point(cls, spins):
        """
        点クラスターの相関関数
        """
        return sum(spins)/len(spins)

    @classmethod
    def _pair(cls, spins):
        """
        第一隣接の二体クラスターの相関関数
        """
        sp_gp = cls._nn_combi(spins)
        #pair_gp = itertools.combinations(sp_gp, 2)
        pair_gp = MakePattern.nCrList(sp_gp, 2)
        pair = [y for x in pair_gp for y in MakePattern.make_tree(*x)]
        return cls.sum_spins(pair)

    @classmethod
    def _pair_nn(cls, spins):
        """
        第二隣接の二体クラスターの相関関数
        """
        sp_gp = cls._nn_combi(spins)
        return cls.sum_spins(sp_gp)

    @classmethod
    def _tri(cls, spins):
        """
        第一隣接のみの正三角形クラスターの相関関数
        """
        sp_gp = cls._nn_combi(spins)
        return cls.sum_spins(MakePattern.make_tree(*sp_gp))

    @classmethod
    def _tri_nn(cls, spins):
        """
        第二隣接を含む二等辺三角形クラスターの相関関数
        """
        sp_gp = cls._nn_combi(spins)
        quad = [x[0] + x[1] for x in MakePattern.nCrList(sp_gp, 2)]
        tri = [y for x in quad for y in MakePattern.nCrList(x, 3)]
        return cls.sum_spins(tri)

    @classmethod
    def _tetra(cls, spins):
        """
        四面体クラスターの相関関数
        """
        return cls._pair(spins) * cls._hex(spins)

    @classmethod
    def _squa(cls, spins):
        """
        正方形クラスターの相関関数
        """
        return cls._pair_nn(spins) * cls._hex(spins)

    @classmethod
    def _penta(cls, spins):
        """
        五体クラスターの相関関数
        """
        return cls._point(spins) * cls._hex(spins)

    @classmethod
    def _hex(cls, spins):
        """
        六体クラスターの相関関数
        """
        return cls.times_spins(spins)

    @classmethod
    def print_matrix(cls):
        """
        AAAA, AAAB, AABB, ABBB, BBBBの5つを縦軸に取った相関関数xiの行列を表示
        """
        HEX = [1, 1, 1, 1, 1, 1]
        PENTA = [-1, 1, 1, 1, 1, 1]
        QUAD = [-1, 1, -1, 1, 1, 1]  # 立体
        SQUA = [-1, -1, 1, 1, 1, 1]  # 平面
        TRI = [-1, 1, -1, 1, -1, 1] # 第一
        TRI_NN = [-1, -1, -1, 1, 1, 1] # 第二
        PAIR = [-1, -1, -1, 1, -1, 1] # 第一
        PAIR_NN = [-1, -1, -1, -1, 1, 1] # 第二
        POINT = [-1, -1, -1, -1, -1, 1]
        NULL_C = [-1, -1, -1, -1, -1, -1]
        spins_list = [HEX, PENTA, QUAD, SQUA, TRI,
                      TRI_NN, PAIR, PAIR_NN, POINT, NULL_C]
        clusters = [cls._null, cls._point, cls._pair, cls._pair_nn, cls._tri,
                    cls._tri_nn, cls._tetra, cls._squa, cls._penta, cls._hex]
        print("Octahedron")
        for spins in spins_list:
            xi = ["{0:.3f}".format(c(spins)) for c in clusters]
            print(" ".join(xi))


class Tetrahedron(object):
    """
    ECIを引数にobjectを生成
    tetraのspin配列を代入することでenergyを算出できる
    """
    def __init__(self, ecis):
        self.ecis = ecis


    @staticmethod
    def _null(_):
        """
        nullクラスター
        """
        return 1

    @classmethod
    def _point(cls, spins):
        """
        点クラスターの相関関数
        """
        return spins.mean()

    @classmethod
    def _pair(cls, spins):
        """
        第一隣接の二体クラスターの相関関数
        """
        idx = np.fromiter(chain.from_iterable(
            combinations(range(len(spins)), 2)), np.int).reshape(-1, 2)
        return (spins[idx[:, 0]] * spins[idx[:, 1]]).mean()

    @classmethod
    def _tri(cls, spins):
        """
        第一隣接のみの正三角形クラスターの相関関数
        """
        return cls._tetra(spins) * cls._point(spins)

    @classmethod
    def _tetra(cls, spins):
        """
        四面体クラスターの相関関数
        """
        return spins.prod()

    @classmethod
    def all(cls, spins):
        """
        [null, point, pair, tri, tetra]の相関関数をlistにしてretrun
        """
        clusters = [cls._null, cls._point, cls._pair, cls._tri, cls._tetra]
        return np.array([c(spins) for c in clusters])

    @classmethod
    def average(cls, spins_array):
        """
        spins_arrayから相関関数を全て計算してその平均をreturn
        """
        return np.array([cls.all(s) for s in spins_array]).mean(0)
    @classmethod
    def print_matrix(cls):
        """
        AAAA, AAAB, AABB, ABBB, BBBBの5つを縦軸に取った相関関数xiの行列を表示
        """
        AAAA = [1, 1, 1, 1]
        AAAB = [1, 1, 1, -1]
        AABB = [1, 1, -1, -1]
        ABBB = [1, -1, -1, -1]
        BBBB = [-1, -1, -1, -1]
        spins_array = np.array([AAAA, AAAB, AABB, ABBB, BBBB])
        clusters = [cls._null, cls._point, cls._pair, cls._tri, cls._tetra]

        print("Tetrahedron")
        for spins in spins_array:
            xi = ["{0:.3f}".format(c(spins)) for c in clusters]
            print(" ".join(xi))

if __name__ == '__main__':
    #for i in range(3000):
    #    Tetrahedron.print_matrix()
    Octahedron.print_matrix()
    Tetrahedron.print_matrix()

