#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Cluster密度を計算
ToDo:iterationで最安定に到達するかをチェックする
"""
from __future__ import division
import pickle
import random
import numpy as np
from itertools import product, chain, combinations

def main():
    """
    execute test
    """
    print('just test')
    cell = CellFCC(10, arrange='L10')
    print(cell.from_origin(np.array([[0, 0, 0, 0]]), 1))


class MonteCarlo(object):
    """
    ECIs, cell, 温度Tをセットに読み込んで, MonteCarloのiteration計算を行う
    """
    def __init__(self, ecis, cell, T):
        self.cell = cell
        self.kb = 8.6171e-5  # eV/K
        self.T = T
        self.beta = self.kb * self.T
        self.ecis = ecis

    def delta_e(self):
        """
        spin flipによる各サイト毎のエネルギーの変化量
        """
        delta_xi_tetra = self.cell.get_delta_xi_tetra()
        return (delta_xi_tetra * self.ecis).sum(axis=1)

    def _delta_e(self, site_id):
        delta_xi_tetra = self.cell._get_delta_xi_tetra(site_id)
        return (delta_xi_tetra * self.ecis).sum(axis=1)

    def _delta_e_TO(self, site_id):
        delta_xi_TO = self.cell._get_delta_xi_TO(site_id)
        return (delta_xi_TO * self.ecis).sum(axis=1)

    def transition_probability(self, de_array):
        """
        エネルギー差の行列から遷移確率を求める
        文献のtype 2
        """
        # 以下は温度を無視した判定処理
        # ToDo:分離して残す
        # de_array[de_array >= 0] = 0
        # de_array[de_array < 0] = 1
        # return de_array

        w = np.exp(-de_array/self.beta)
        return w/(1+w)

    def iteration_lump(self, steps):
        """
        site毎にわけずに行列処理で一挙にflipの判定を行う
        この方法では単振動してしまう場合がある
        """
        print("initial energy")
        print((self.cell.get_tetra_xi_glob() * self.ecis).sum())
        for _ in range(steps):
            de = self.delta_e()
            p = self.transition_probability(de)
            rand = np.random.rand(p.size)
            judge = rand < p
            judge = judge.reshape(4, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 4)
            self.cell.cell[judge] *= -1
        print((self.cell.get_tetra_xi_glob() * self.ecis).sum())
        return self.cell.cell

    def iteration(self, steps):
        """
        site0, 1, 2, 3の四種を別けてflip判定
        """
        print("initial energy")
        print((self.cell.get_tetra_xi_glob() * self.ecis).sum())
        for _ in range(steps):
            for j in range(4):
                de = self.delta_e()
                p = self.transition_probability(de)
                rand = np.random.rand(p.size)
                judge = rand < p

                judge = judge.reshape(4, -1).T.reshape(
                    self.cell.size, self.cell.size, self.cell.size, 4)
                judge2 = np.array([False]*self.cell.size**3*4).reshape(
                    self.cell.size, self.cell.size, self.cell.size, 4)
                judge2[:, :, :, j] = judge[:, :, :, j]
                self.cell.cell[judge2] *= -1
                conc_ene = [self.cell.get_conc(1),
                            (self.cell.get_tetra_xi_glob() * self.ecis).sum()]
                print("{0[0]} {0[1]}".format(conc_ene))
        return self.cell.cell

    def _iteration(self, step):
        """
        """
        #print("initial energy")
        #print((self.cell.get_tetra_xi_glob() * self.ecis).sum())
        conc_ene = [[self.cell.get_conc(1),
                     (self.cell.get_tetra_xi_glob() * self.ecis).sum()]]
        for _ in range(step):
            for j in range(4):
                de = self._delta_e(j)
                p = self.transition_probability(de)
                rand = np.random.rand(p.size)
                judge = (rand < p).reshape(
                    self.cell.size, self.cell.size, self.cell.size)
                judge2 = np.array([False]*self.cell.size**3*4).reshape(
                    self.cell.size, self.cell.size, self.cell.size, 4)
                judge2[:, :, :, j] = judge
                self.cell.cell[judge2] *= -1
                conc_ene.append(
                    [self.cell.get_conc(1),
                     (self.cell.get_tetra_xi_glob() * self.ecis).sum()])
        return conc_ene

    def _iterationTO(self, step):
        """
        TOによる
        """
        conc_ene = [[self.cell.get_conc(1),
                     (self.cell.get_xi_TO_glob() * self.ecis).sum()]]
        print("initial")
        print(conc_ene[-1])
        print("start")
        for _ in range(step):
            for j in range(4):
                de = self._delta_e_TO(j)
                p = self.transition_probability(de)
                rand = np.random.rand(p.size)
                judge = (rand < p).reshape(
                    self.cell.size, self.cell.size, self.cell.size)
                judge2 = np.array([False]*self.cell.size**3*4).reshape(
                    self.cell.size, self.cell.size, self.cell.size, 4)
                judge2[:, :, :, j] = judge
                self.cell.cell[judge2] *= -1
                conc_ene.append(
                    [self.cell.get_conc(1),
                     (self.cell.get_xi_TO_glob() * self.ecis).sum()])
                print(conc_ene[-1])
        return conc_ene

    def _iterationTO_reserved_atoms(self, step):
        """
        TOによる
        粒子数を保存
        """
        conc_ene = [[self.cell.get_conc(1),
                     (self.cell.get_xi_TO_glob() * self.ecis).sum()]]
        print("initial")
        print(conc_ene[-1])
        print("start")
        for _ in range(step):
            for j in range(4):
                de = self._delta_e_TO(j)
                p = self.transition_probability(de)
                rand = np.random.rand(p.size)
                judge = (rand < p).reshape(
                    self.cell.size, self.cell.size, self.cell.size)
                judge2 = np.array([False]*self.cell.size**3*4).reshape(
                    self.cell.size, self.cell.size, self.cell.size, 4)
                judge2[:, :, :, j] = judge
                m2p = judge2 * (self.cell.cell == -1)
                p2m = judge2 * (self.cell.cell == 1)
                is_m2p_major = m2p.sum() > p2m.sum()
                major = {True: m2p, False: p2m}[is_m2p_major]
                minor = {True: p2m, False: m2p}[is_m2p_major]
                diff = major.sum() - minor.sum()
                # print(diff)
                if diff != 0:
                    fix = np.random.choice(
                        range(major.sum()), diff, replace=False)
                    major[tuple(np.array(np.where(major))[:, fix])] = False

                self.cell.cell[major] *= -1
                self.cell.cell[minor] *= -1

                conc_ene.append(
                    [self.cell.get_conc(1),
                     (self.cell.get_xi_TO_glob() * self.ecis).sum()])
                print(conc_ene[-1])
        return conc_ene

class Octahedron(object):
    """
    四面体と異なってサイトの位置関係で相関関数が異なる
    spinsは(0,1) (2,3) (4,5)が第二隣接の組み合わせとして入力する
    spinsを第二隣接同士の各組み合わせで平均をとって[-1/0/1]*3のlistに
    変換して処理した方が早いかもしれない
    パターンとして3^3=27通りを分類分けすることになる
    """
    def __init__(self):
        pass

    @staticmethod
    def _null(spins):
        """
        nullクラスター
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return np.array([1]*spins.shape[0])

    @staticmethod
    def _point(spins):
        """
        点クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins.mean(2).mean(1)

    @staticmethod
    def _pair(spins):
        """
        第一隣接の二体クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        # fancy indexを使う
        # 第二隣接以外の4体のspinから1体を選ぶ全ての組み合わせ
        s = spins.reshape(-1, 6, 3, 2)

        # 3C2
        idx1 = np.fromiter(chain.from_iterable(
            combinations(range(3), 2)), np.int).reshape(-1, 2)
        # 2*2
        idx2 = np.fromiter(chain.from_iterable(
            product(range(2), repeat=2)), np.int).reshape(-1, 2)

        return (s[:, :, idx1[:, 0]][:, :, :, idx2[:, 0]] *
                s[:, :, idx1[:, 1]][:, :, :, idx2[:, 1]]).mean(
                    axis=3).mean(axis=2).mean(axis=1)

    @classmethod
    def _pair_nn(cls, spins):
        """
        第二隣接の二体クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins.reshape(-1, 6, 3, 2).prod(axis=3).mean(axis=2).mean(axis=1)

    @classmethod
    def _tri(cls, spins):
        """
        第一隣接のみの正三角形クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        s = spins.reshape(-1, 6, 3, 2)
        # 2^3
        idx = np.fromiter(chain.from_iterable(
            product(range(2), repeat=3)), np.int).reshape(-1, 3)
        return s[:, :, [0, 1, 2], idx].reshape(
            -1, 6, 8, 3).prod(axis=3).mean(axis=2).mean(axis=1)

    @classmethod
    def _tri_nn(cls, spins):
        """
        第二隣接を含む二等辺三角形クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        s = spins.reshape(-1, 6, 3, 2)
        idx = np.fromiter(chain.from_iterable(
            combinations(range(3), r=2)), np.int).reshape(-1, 2)
        s = s[:, :, idx, :].reshape(-1, 6, 3, 4)
        idx = np.fromiter(chain.from_iterable(
            combinations(range(4), r=3)), np.int).reshape(-1, 3)
        return s[:, :, :, idx].reshape(
            -1, 6, 12, 3).prod(axis=3).mean(axis=2).mean(axis=1)

    @classmethod
    def _tetra(cls, spins):
        """
        四面体クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        s = spins.reshape(-1, 6, 3, 2)

        # 3C2
        idx1 = np.fromiter(chain.from_iterable(
            combinations(range(3), 2)), np.int).reshape(-1, 2)
        # 2*2
        idx2 = np.fromiter(chain.from_iterable(
            product(range(2), repeat=2)), np.int).reshape(-1, 2)
        return ((s[:, :, idx1[:, 0]][:, :, :, idx2[:, 0]] *
                 s[:, :, idx1[:, 1]][:, :, :, idx2[:, 1]]).mean(
                    axis=3).mean(axis=2) * spins.prod(axis=2)).mean(axis=1)
        return cls._pair(spins) * cls._hex(spins)


    @classmethod
    def _squa(cls, spins):
        """
        正方形クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """

        return (spins.reshape(-1, 6, 3, 2).prod(axis=3).mean(axis=2) *
                spins.prod(axis=2)).mean(axis=1)
        return cls._pair_nn(spins) * cls._hex(spins)

    @classmethod
    def _penta(cls, spins):
        """
        五体クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return cls._point(spins) * cls._hex(spins)

    @classmethod
    def _hex(cls, spins):
        """
        六体クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins.prod(axis=2).mean(axis=1)

    @classmethod
    def get_averaged_xi(cls, spins):
        """
        10種の相関関数の平均値をarrayにしてretrun
        """
        return cls.get_xi(spins).mean(axis=0)

    @classmethod
    def get_xi(cls, spins):
        """
        10種の相関関数をarrayにしてretrun
        attributes
            spins: N×8×4
        return
            xi: N×10個のxi (各サイト毎の値を出力)
        """
        clusters = [cls._null, cls._point, cls._pair, cls._pair_nn, cls._tri,
                    cls._tri_nn, cls._tetra, cls._squa, cls._penta, cls._hex]
        return np.r_[[c(spins) for c in clusters]].T

    @classmethod
    def print_matrix(cls):
        """
        AAAA, AAAB, AABB, ABBB, BBBBの5つを縦軸に取った相関関数xiの行列を表示
        """
        AAAAAA = np.array([1, 1, 1, 1, 1, 1] * 6).reshape(1, 6, 6)
        AAAAAB = np.array([1, 1, 1, 1, 1, -1] * 6).reshape(1, 6, 6)
        AAAABB = np.array([1, 1, 1, 1, -1, -1] * 6).reshape(1, 6, 6)
        AAABAB = np.array([1, 1, 1, -1, 1, -1] * 6).reshape(1, 6, 6)
        AAABBB = np.array([1, 1, 1, -1, -1, -1] * 6).reshape(1, 6, 6)
        ABABAB = np.array([1, -1, 1, -1, 1, -1] * 6).reshape(1, 6, 6)
        AABBBB = np.array([1, 1, -1, -1, -1, -1] * 6).reshape(1, 6, 6)
        ABABBB = np.array([1, -1, 1, -1, -1, -1] * 6).reshape(1, 6, 6)
        ABBBBB = np.array([1, -1, -1, -1, -1, -1] * 6).reshape(1, 6, 6)
        BBBBBB = np.array([-1, -1, -1, -1, -1, -1] * 6).reshape(1, 6, 6)
        spins_array = np.array([AAAAAA, AAAAAB, AAABAB, AAAABB, ABABAB,
                                AAABBB, ABABBB, AABBBB, ABBBBB, BBBBBB])

        print("Octahedron")
        for spins in spins_array:
            xi = ["{0:.3f}".format(x) for x in cls.get_averaged_xi(spins)]
            print(" ".join(xi))


class Tetrahedron(object):
    """
    ECIを引数にobjectを生成
    四体のspin配列のarray(n*4行列)を代入することでenergy及び相関関数を算出
    FCCでは一つの格子点中心に8つの四面体が存在する
    N*8*4のmatrixでスピン情報を受け取る
    spin flipによる相関関数の変化量が計算ができるよう各点毎の相関関数を
    returnする
    """
    CuAu = [[1.77281, -3693.31368, 1932616.12707],
            [0.35925, -1576.20925, 1254312.59500],
            [-0.00413, -104.24365, 202856.45586],
            [0.01325, -36.57425, 27836.99500],
            [0.01381, -30.42018, 16596.62707]]

    def __init__(self, conc):
        self.conc = conc
        self.a = 6.7318 + 0.9776 * conc

    @property
    def ecis(self):
        """
        eciをreturn
        """
        eci = (self.CuAu *
               np.array([1, self.a**(-3.5), self.a**(-7)])).sum(axis=1)
        eci[1] = -0.011 * 8
        Ry2eV = 13.6058
        eci = eci * Ry2eV
        return eci
        #return (self.CuAu *
        #        np.array([1, self.a**(-3.5), self.a**(-7)])).sum(axis=1)

    @staticmethod
    def _null(spins):
        """
        nullクラスター
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        return np.array([1]*spins.shape[0])

    @staticmethod
    def _point(spins):
        """
        点クラスターの相関関数
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        return spins.mean(2).mean(1)

    @staticmethod
    def _pair(spins):
        """
        第一隣接の二体クラスターの相関関数
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        # fancy indexを使う
        # 4体のspinから2体を選ぶ全ての組み合わせ 4C2
        idx = np.fromiter(chain.from_iterable(
            combinations(range(4), 2)), np.int).reshape(-1, 2)
        # 各組み合わせで積を取って平均をとる
        return spins[:, :, idx].prod(axis=3).mean(axis=2).mean(axis=1)

    @staticmethod
    def _tri(spins):
        """
        第一隣接正三角形クラスターの相関関数
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        # tri = point * tetraによって算出できる
        #return (spins.mean(axis=2) * spins.prod(axis=2)).mean(axis=1)
        # しかし、spin = 1/0に取り直すことができなくなるのでちゃんと計算する
        # 4体のspinから3体を選ぶ全ての組み合わせ 3C2
        idx = np.fromiter(chain.from_iterable(
            combinations(range(4), r=3)), np.int).reshape(-1, 3)
        return spins[:, :, idx].prod(axis=3).mean(axis=2).mean(axis=1)

    @staticmethod
    def _tetra(spins):
        """
        四面体クラスターの相関関数
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        return spins.prod(axis=2).mean(axis=1)


    @classmethod
    def get_averaged_xi(cls, spins):
        """
        [null, point, pair, tri, tetra]の相関関数の平均値をarrayにしてretrun
        """
        return cls.get_xi(spins).mean(axis=0)

    @classmethod
    def get_xi(cls, spins):
        """
        [null, point, pair, tri, tetra]の相関関数をarrayにしてretrun
        attributes
            spins: N×8×4
        return
            xi: N×5個のxi (各サイト毎の値を出力)
        """
        clusters = [cls._null, cls._point, cls._pair, cls._tri, cls._tetra]
        return np.r_[[c(spins) for c in clusters]].T

    @classmethod
    def print_matrix(cls):
        """
        AAAA, AAAB, AABB, ABBB, BBBBの5つを縦軸に取った相関関数xiの行列を表示
        """
        AAAA = np.array([1, 1, 1, 1] * 8).reshape(1, 8, 4)
        AAAB = np.array([1, 1, 1, -1] * 8).reshape(1, 8, 4)
        AABB = np.array([1, 1, -1, -1] * 8).reshape(1, 8, 4)
        ABBB = np.array([1, -1, -1, -1] * 8).reshape(1, 8, 4)
        BBBB = np.array([-1, -1, -1, -1] * 8).reshape(1, 8, 4)
        spins_array = np.array([AAAA, AAAB, AABB, ABBB, BBBB])

        print("Tetrahedron")
        for spins in spins_array:
            xi = ["{0:.3f}".format(x) for x in cls.get_averaged_xi(spins)]
            print(" ".join(xi))

    @classmethod
    def print_energies(cls, ecis):
        """
        ecisに基づいてAAAA, AAAB, AABB, ABBB, BBBBの5種のエネルギーをprint
        """
        AAAA = np.array([1, 1, 1, 1] * 8).reshape(1, 8, 4)
        AAAB = np.array([1, 1, 1, -1] * 8).reshape(1, 8, 4)
        AABB = np.array([1, 1, -1, -1] * 8).reshape(1, 8, 4)
        ABBB = np.array([1, -1, -1, -1] * 8).reshape(1, 8, 4)
        BBBB = np.array([-1, -1, -1, -1] * 8).reshape(1, 8, 4)
        spins_array = np.array([AAAA, AAAB, AABB, ABBB, BBBB])

        print("Energies")
        for spins in spins_array:
            energy = (cls.get_averaged_xi(spins) * ecis).sum()
            print(energy)


class QuadSite(object):
    """
    4点のサイトの情報を記述する
    spins
    sites(不使用、memo)
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


class CellFCC(object):
    """
    np.array object of 3dim×4sites
    size: cellの大きさ (size^3)
    arrange: 'random', 'L12', 'L10'を指定
    conc: A, B 二元系におけるAの濃度の初期値 (spinはA=1, B=-1とする)
    """
    # 原点からみた最隣接原子の座標
    NEAREST = [[0, 0, 0, 1], [0, 0, 0, 2], [0, -1, 0, 1], [-1, 0, 0, 2],
               [0, 0, 0, 3], [0, -1, 0, 3], [-1, -1, 0, 3], [-1, 0, 0, 3],
               [0, 0, -1, 1], [0, 0, -1, 2], [0, -1, -1, 1], [-1, 0, -1, 2]]
    # memo 上記の実座標 簡単のため2倍して表記
    COORDINATE_NEAREST = [[0, 1, 1], [1, 0, 1], [0, -1, 1], [-1, 0, 1],
                          [1, 1, 0], [1, -1, 0], [-1, -1, 0], [-1, 1, 0],
                          [0, 1, -1], [1, 0, -1], [0, -1, -1], [-1, 0, -1]]
    # NEARESTを原点＋3つ選んで作る8つの四面体 (原点は省略)
    TETRA_POS = [[0, 1, 4], [3, 0, 7], [2, 3, 6], [1, 2, 5],
                 [8, 9, 4], [11, 8, 7], [10, 11, 6], [9, 10, 5]]

    # 原点から見た第二隣接原子の座標
    NEXT_N = [[1, 0, 0, 0], [-1, 0, 0, 0],
              [0, 1, 0, 0], [0, -1, 0, 0],
              [0, 0, 1, 0], [0, 0, -1, 0]]

    # 第二隣接一つと4つの最隣接からなる八面体
    # その中の NEAREST 4つの組み合わせ
    # 原点と原点に対しての第二隣接は省略 (省略の第二隣接は NEXT_N の順番と対応)
    # また OCTA_POS 中の[前2、後2]も第二隣接の関係にある (例: 1,9 及び 4,5)
    OCTA_POS = [[1, 9, 4, 5], [3, 11, 6, 7],
                [0, 8, 4, 7], [2, 10, 5, 6],
                [0, 2, 1, 3], [8, 10, 9, 11]]

    def __init__(self, size, arrange='random', conc=0.5):
        self.size = size
        mode = {'random': QuadSite.random, 'L10': QuadSite.L10,
                'L12_A': QuadSite.L12_A, 'L12_B': QuadSite.L12_B,}[arrange]
        self.cell = np.array([[[mode(conc).spins for z in range(size)]
                               for y in range(size)] for x in range(size)])

        self.TETRA = [[[0, 0, 0, 0]] + [self.NEAREST[y] for y in x]
                      for x in self.TETRA_POS]
        self.OCTA = [[[0, 0, 0, 0]] + [nn] + [self.NEAREST[y] for y in x]
                     for x, nn in zip(self.OCTA_POS, self.NEXT_N)]

        # 原点を変更したとき、4つのsiteから見たquad_siteの相対位置
        # 上から順に原点に移すsite_id=0,1,2,3及び
        # 右から順に新しい原点から見た相対的なsite_id=0,1,2,3を示している
        # 例えばsite_id=1を原点としてsite_id=2を見た場合、
        # 座標を+[0,0,1,1]移動してspinを参照する
        self.TRANS = np.array(
            [[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
             [[0, 0, 0, 1], [0, 1, 1, -1], [0, 0, 1, 1], [0, 1, 0, -1]],
             [[0, 0, 0, 2], [0, 0, 1, 2], [1, 0, 1, -2], [1, 0, 0, -2]],
             [[0, 0, 0, 3], [0, 1, 0, 1], [1, 0, 0, -1], [1, 1, 0, -3]]])

        # 3次元のすべての座標
        self.coords = list(product(range(0, self.size), repeat=3))
        # site idは0のみ
        self.coords_id0 = [list(x) + [0] for x in self.coords]

    def make_poscar(self, fname):
        """
        POSCARを作成する
        """
        poss1 = self.get_real_coordinate(
            np.argwhere(self.cell == 1)) / self.size
        poss2 = self.get_real_coordinate(
            np.argwhere(self.cell != 1)) / self.size
        num1 = poss1.shape[0]
        num2 = poss2.shape[0]
        lines = "mc\n"
        lines += "{0}\n".format(self.size*2.5)
        lines += "1.00 0 0\n0 1.00 0\n0 0 1.00\n"
        lines += "Al Cu\n"
        lines += "{0} {1}\n".format(num1, num2)
        lines += "Direct\n"
        for site in poss1:
            lines += " ".join([str(x) for x in site]) + "\n"
        for site in poss2:
            lines += " ".join([str(x) for x in site]) + "\n"
        with open(fname, 'w') as wfile:
            wfile.write(lines)

    def make_poscar_order(self, fname, ecis):
        """
        POSCARを作成する
        """
        delta_xi_TO = np.hstack([self._get_delta_xi_TO(x)
                                 for x in range(4)]).reshape(-1, 4, 11)
        de = (delta_xi_TO * ecis).sum(axis=2).reshape(
            self.size, self.size, self.size, 4)

        p1 = (self.cell == 1) * (de > 0.032)
        p2 = (self.cell != 1) * (de > 0.032)
        p3 = de <= 0.32
        poss1 = self.get_real_coordinate(
            np.argwhere(p1)) / self.size
        poss2 = self.get_real_coordinate(
            np.argwhere(p2)) / self.size
        poss3 = self.get_real_coordinate(
            np.argwhere(p3)) / self.size
        num1 = poss1.shape[0]
        num2 = poss2.shape[0]
        num3 = poss3.shape[0]
        lines = "mc\n"
        lines += "{0}\n".format(self.size*2.5)
        lines += "1.00 0 0\n0 1.00 0\n0 0 1.00\n"
        lines += "Al Cu Ni\n"
        lines += "{0} {1} {2}\n".format(num1, num2, num3)
        lines += "Direct\n"
        for site in poss1:
            lines += " ".join([str(x) for x in site]) + "\n"
        for site in poss2:
            lines += " ".join([str(x) for x in site]) + "\n"
        for site in poss3:
            lines += " ".join([str(x) for x in site]) + "\n"
        with open(fname, 'w') as wfile:
            wfile.write(lines)

    def get_real_coordinate(self, coords):
        """
        get coordinate in real space
        """
        trans = np.array([[0, 0, 0], [0, 0.5, 0.5],
                          [0.5, 0, 0.5], [0.5, 0.5, 0]])
        r_coords = coords[:, :3] + trans[coords[:, 3]]
        return r_coords


    def save_cell(self, title):
        fname = title + ".pickle"
        with open(fname, 'wb') as wbfile:
            pickle.dump(self, wbfile)

    @staticmethod
    def load_cell(pickle_path):
        with open(pickle_path, 'rb') as rbfile:
            cell = pickle.load(rbfile)
            return cell

    def get_orderparam(self, size):
        """
        order parameterをreturnする
        order parameterの評価方法は澤田さんのメモを利用
        """
        latt = np.fromiter(chain.from_iterable(
            product(range(2*size+1), repeat=3)), np.int).reshape(-1, 3)
        idx = latt.sum(axis=1) % 2 == 0
        latt0 = np.c_[latt[idx], [0]*latt[idx].shape[0]]
        latt1 = latt0 + np.array([0, 0, 0, 1])
        latt2 = latt0 + np.array([0, 0, 0, 2])
        latt3 = latt0 + np.array([0, 0, 0, 3])
        latt4 = latt0 + np.array([1, 0, 0, 0])
        latt5 = latt4 + np.array([0, 0, 0, 1])
        latt6 = latt4 + np.array([0, 0, 0, 2])
        latt7 = latt4 + np.array([0, 0, 0, 3])

        p0 = self.cell[tuple(latt0.T)].mean()
        p1 = self.cell[tuple(latt1.T)].mean()
        p2 = self.cell[tuple(latt2.T)].mean()
        p3 = self.cell[tuple(latt3.T)].mean()
        p4 = self.cell[tuple(latt4.T)].mean()
        p5 = self.cell[tuple(latt5.T)].mean()
        p6 = self.cell[tuple(latt6.T)].mean()
        p7 = self.cell[tuple(latt7.T)].mean()

        s1 = (p0 + p1 - p2 - p3 + p4 + p5 - p6 - p7) / 8
        s2 = (p0 - p1 + p2 - p3 + p4 - p5 - p6 + p7) / 8
        s3 = (p0 - p1 - p2 + p3 + p4 - p5 - p6 + p7) / 8
        #print(p0, p1, p2, p3, p4, p5, p6, p7)
        return s1, s2, s3

    def get_xi_TO(self, site_id):
        """
        tetraとoctaの2種を合わせた各サイト毎のxiをreturn
        [null, point, pair, pair_nn, tri, tri_nn, tetra, tetra_octa,
         squa, penta, hex]の11種
        attributes
            spins: N×8×4
        return
            xi: N×5個のxi (各サイト毎の値を出力)
        """
        pass



    def get_xi_TO_glob(self):
        """
        tetra, octaを合わせて出力
        """
        octa = self.get_octa_xi_glob()
        tetra = self.get_tetra_xi_glob()
        #if not all((octa[[0, 1, 2, 4]] - tetra[:4]) ** 2 < 1e-8 ):
        #    print("xi of tetra and octa is different")
        return np.array([tetra[0], tetra[1], tetra[2], octa[3], tetra[3],
                         octa[5], tetra[4], octa[6], octa[7], octa[8], octa[9]])

    def __get_octa_xi(self, coord):
        """
        不要
        座標(coord)周りの6つのoctahedronクラスターの相関関数を計算して
        その平均をreturnする
        attributes
            coord: (x, y, z, site_id)
        """
        # coord 周りの octahedron の座標の array (6組*6点*座標4)
        octa = (np.array(coord[:3]+[0]) + self.OCTA).reshape(36, 4)

        # 八面体クラスターのspin(6つ)の組み合わせ*6組
        spins = self.from_origin(octa, coord[3])
        spins = spins.reshape(1, 6, 6)
        return Octahedron.get_averaged_xi(spins)

    def get_octa_xi_glob(self):
        """
        系内全域の相関関数をreturn
        """
        octa = (np.array(self.coords_id0).reshape(self.size**3, 1, 1, 4) +
                self.OCTA).reshape(self.size**3*6*6, 4)

        # site id 4種 行列として結合
        spins = np.r_[[self.from_origin(octa, i) for i in range(4)]]
        spins = spins.reshape(self.size**3*4, 6, 6)
        return Octahedron.get_averaged_xi(spins)

    def __get_tetra_xi(self, coord):
        """
        不要
        座標(coord)周りの8つのtetragonalクラスターの相関関数を計算して
        その平均をreturnする
        attributes
            coord: (x, y, z, site_id)
        """
        # coord周りのtetrahedronの座標のarray (8組*4点*座標4)
        tetra = (np.array(coord[:3]+[0]) + self.TETRA).reshape(32, 4)

        # 四面体クラスターのspin(4つ)の組み合わせ*8組
        spins = self.from_origin(tetra, coord[3])
        spins = spins.reshape(1, 8, 4)
        return Tetrahedron.get_averaged_xi(spins)

    def get_tetra_xi_glob(self):
        """
        系内全域の相関関数をreturn
        """
        tetra = (np.array(self.coords_id0).reshape(self.size**3, 1, 1, 4) +
                 self.TETRA).reshape(self.size**3*8*4, 4)

        # site id 4種 行列として結合
        spins = np.r_[[self.from_origin(tetra, i) for i in range(4)]]
        spins = spins.reshape(self.size**3*4, 8, 4)
        return Tetrahedron.get_averaged_xi(spins)

    def get_delta_xi_tetra(self):
        """
        spinがフリップすることによる相関関数の変化量を各サイト毎、
        すべてのサイトに渡って算出する
        8つのtetrahedronのみを計算
        return
            N * 5
        """
        tetra = (np.array(self.coords_id0).reshape(self.size**3, 1, 1, 4) +
                 self.TETRA).reshape(self.size**3*8*4, 4)

        # site id 4種を行列として結合
        spins = np.r_[[self.from_origin(tetra, i) for i in range(4)]]
        spins = spins.reshape(self.size**3*4, 8, 4)

        flipped = np.copy(spins)
        flipped[:, :, 0] = spins[:, :, 0] * -1
        return (Tetrahedron.get_xi(flipped) -
                Tetrahedron.get_xi(spins))

    def _get_delta_xi_tetra(self, site_id):
        """
        spinがフリップすることによる相関関数の各座標毎の変化量を
        すべての座標に渡って算出する
        8つのtetrahedronのみを計算
        return
            N * 5
        """
        tetra = (np.array(self.coords_id0).reshape(self.size**3, 1, 1, 4) +
                 self.TETRA).reshape(self.size**3*8*4, 4)
        spins = self.from_origin(tetra, site_id).reshape(self.size**3, 8, 4)
        flipped = np.copy(spins)
        flipped[:, :, 0] *= -1
        return Tetrahedron.get_xi(flipped) - Tetrahedron.get_xi(spins)

    def _get_delta_xi_TO(self, site_id):
        """
        spinがフリップすることによる相関関数の各座標毎の変化量を
        すべての座標に渡って算出する
        tetrahedron and octahederon
        return
            N * 5
        """
        tetra = (np.array(self.coords_id0).reshape(self.size**3, 1, 1, 4) +
                 self.TETRA).reshape(self.size**3*8*4, 4)
        spins = self.from_origin(tetra, site_id).reshape(self.size**3, 8, 4)
        flipped = np.copy(spins)
        flipped[:, :, 0] *= -1
        dxi_tet = Tetrahedron.get_xi(flipped) - Tetrahedron.get_xi(spins)

        octa = (np.array(self.coords_id0).reshape(self.size**3, 1, 1, 4) +
                self.OCTA).reshape(self.size**3*6*6, 4)
        spins = self.from_origin(octa, site_id).reshape(self.size**3, 6, 6)
        flipped = np.copy(spins)
        flipped[:, :, 0] *= -1
        dxi_octa = Octahedron.get_xi(flipped) - Octahedron.get_xi(spins)
        return np.r_[[dxi_tet[:, 0], dxi_tet[:, 1], dxi_tet[:, 2], dxi_octa[:, 3],
                      dxi_tet[:, 3], dxi_octa[:, 5], dxi_tet[:, 4], dxi_octa[:, 6],
                      dxi_octa[:, 7], dxi_octa[:, 8], dxi_octa[:, 9]]].T



    def _dxi_tetra(self,coords):
        """
        """
        pass

    def spin_flip_dxi(self, coord):
        """
        tetraとoctaの2種をreturn
        """
        return self.get_delta_xi_tetra()

    def from_origin(self, coords, site_id=0):
        """
        任意のsite_id中心から見た相対座標のlist(coords)を引数にして、
        その座標におけるspinのリストをreturn
        端でsize+2となる点が現れるのでそれを補正する為に%= self.sizeする
        attributes
            coords: site_idからの相対座標
            site_id: 原点に据えるQuadのsite id
        """
        rpos = coords + self.TRANS[site_id][coords[:, 3]]
        rpos[:, 0:3] %= self.size
        return self[tuple(np.transpose(rpos))]

    def get_conc(self, i):
        """
        A元素の格子中の濃度を算出
        """
        spin = [s for x in self.cell for y in x for z in y for s in z]
        return spin.count(i) / len(spin)

    def __getitem__(self, *pos):
        return self.cell.__getitem__(*pos)

if __name__ == '__main__':
    main()


