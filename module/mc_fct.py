#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
monte carlo計算用のモジュール
MonteCarlo: スピンフリップ判定を取り纏める
Octahedron: 八面体クラスターの相関関数 (立方晶表記)
Tetrahedron: 四面体クラスターの相関関数 (立方晶表記)
TetraTetra: 四面体クラスターの相関関数 (正方晶表記)
QuadSite: 4副格子のサイト情報を記述 不要？？
FaceCenterCubic: FCC格子のスピン情報を記憶する
FaceCenterTetragonal: FCT格子のスピン情報を記憶する (a = b ≠ c)
"""
from __future__ import division
import pickle
import random
import numpy as np
from itertools import product, chain, combinations, permutations

def main():
    """
    execute test
    """
    print('just test')
    cell = FaceCenterTetragonal(10, arrange='L10')
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
        delta_xi_tetra = self.cell.get_tetra_dxi_glob()
        return (delta_xi_tetra * self.ecis).sum(axis=1)

    def _delta_e(self, site_id):
        delta_xi_tetra = self.cell.get_tetra_dxi(site_id)
        return (delta_xi_tetra * self.ecis).sum(axis=1)

    def _delta_e_TO(self, site_id):
        delta_xi_TO = self.cell.get_TO_dxi(site_id)
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
                     (self.cell.get_TO_xi_glob() * self.ecis).sum()]]
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
                     (self.cell.get_TO_xi_glob() * self.ecis).sum()])
                print(conc_ene[-1])
        return conc_ene

    def _iterationTO_reserved_atoms(self, step):
        """
        TOによる
        粒子数を保存
        """
        conc_ene = [[self.cell.get_conc(1),
                     (self.cell.get_TO_xi_glob() * self.ecis).sum()]]
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
                     (self.cell.get_TO_xi_glob() * self.ecis).sum()])
                print(conc_ene[-1])
        return conc_ene


class Octahedron(object):
    """
    spinsは(0,1) (2,3)が第三隣接の組み合わせ、
    (4,5)が第四隣接の組み合わせとして入力する
    other variables:
        ### 二体
        _idx_pairR3: R3同士のpairを組にして並べたindex
        _idx_pairR4: R4のpairの組み合わせを示すindex
        _idx_pairR1: R1同士のpairを組にして並べたindex
                     R3 2組の直積
        _idx_pairR2: R2のpairの組み合わせを示すindex
                     R4から一つ、それ以外から一つ選ぶ組み合わせ
        ### 三体
        _idx_triR2R1: R2 二つと R1 一つを辺に持つ三角形
                      R4から一つ、pairR1(R3 2組の直積)から一つを選ぶ組み合わせ
        _idx_triR3R1: R3 一つと R1 二つを辺に持つ三角形
                      (0, 1, 2, 3)から三つ選ぶ
        _idx_triR3R2: R3 一つと R2 二つを辺に持つ三角形
                      ([0, 1], [2, 3]) から一組 (4, 5)から一つ選ぶ
        _idx_triR4R2: R4 二つと R2 二つを辺に持つ三角形
                      (0, 1, 2, 3) から一つと [4, 5]の二つ
        ### 四体 (R*は2つまで記載)
        _idx_tetraR3R1: R1(辺) 四つ R3(対角) 二つ からなる平面四角形
                      六体のspinからpairR4を一つ取り除いたもの
        _idx_tetraR4R3: R2(辺)四つ R4(対角)一つ R3(対角)一つ からなる平面四角形
                          六体のspinからpairR3を一つ取り除いたもの
        _idx_tetraR3R2: R1 二つ R2 三つ R3 一つ からなる四面体
                          六体のspinからpairR2を一つ取り除いたもの
        _idx_tetraR4R2: R1 一つ R2 四つ R4 一つ からなる四面体
                          六体のspinからpairR1を一つ取り除いたもの
        ### 五体 (R*は最大のもののみ記載)
        _idx_pentaR3: R4を含まない五体のクラスター
        _idx_pentaR4: R4を含む五体のクラスター

    """
    _idx_pairR1 = np.fromiter(chain.from_iterable(
        product([0, 1], [2, 3])), np.int).reshape(-1, 2)
    _idx_pairR2 = np.fromiter(chain.from_iterable(
        product([0, 1, 2, 3], [4, 5])), np.int).reshape(-1, 2)
    _idx_pairR3 = np.array([[0, 1], [2, 3]])
    _idx_pairR4 = np.array([[4, 5]])

    _idx_triR2R1 = np.array([list(x) + [y] for x in product([0, 1], [2, 3])
                             for y in [4, 5]]).reshape(-1, 3)
    _idx_triR3R1 = np.fromiter(chain.from_iterable(
        combinations(range(4), 3)), np.int).reshape(-1, 3)
    _idx_triR3R2 = np.array(
        [x + [y] for x in [[0, 1], [2, 3]] for y in [4, 5]]).reshape(-1, 3)
    _idx_triR4R2 = np.array([[x] + [4, 5] for x in range(4)]).reshape(-1, 3)

    _idx_tetraR3R1 = np.array(
        [[x for x in range(6) if x not in y] for y in _idx_pairR4])
    _idx_tetraR4R3 = np.array(
        [[x for x in range(6) if x not in y] for y in _idx_pairR3])
    _idx_tetraR3R2 = np.array(
        [[x for x in range(6) if x not in y] for y in _idx_pairR2])
    _idx_tetraR4R2 = np.array(
        [[x for x in range(6) if x not in y] for y in _idx_pairR1])

    _idx_pentaR3 = np.array(
        [[x for x in range(6) if x != y] for y in [4, 5]])
    _idx_pentaR4 = np.array(
        [[x for x in range(6) if x != y] for y in [0, 1, 2, 3]])

    def __init__(self):
        pass

    @staticmethod
    def _null(spins):
        """
        nullクラスター
        attributes
            spins: N×6×6
        return
            xi: N個のxi (全て1)
        """
        return np.ones(spins.shape[0:-2])

    @staticmethod
    def _point(spins):
        """
        点クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins.mean(-1).mean(-1)

    @classmethod
    def _pairR1(cls, spins):
        """
        第一隣接の二体クラスターの相関関数
        fancy indexを使って算出
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_pairR1].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _pairR2(cls, spins):
        """
        第二隣接の二体クラスターの相関関数
        R2のpairになるよう再配列してprodを取る
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_pairR2].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _pairR3(cls, spins):
        """
        第二隣接の二体クラスターの相関関数
        R2のpairになるよう再配列してprodを取る
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_pairR3].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _pairR4(cls, spins):
        """
        第二隣接の二体クラスターの相関関数
        R2のpairになるよう再配列してprodを取る
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_pairR4].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _triR2R1(cls, spins):
        """
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_triR2R1].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _triR3R1(cls, spins):
        """
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_triR3R1].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _triR3R2(cls, spins):
        """
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_triR3R2].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _triR4R2(cls, spins):
        """
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_triR4R2].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)


    @classmethod
    def _tetraR3R1(cls, spins):
        """
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_tetraR3R1].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _tetraR4R3(cls, spins):
        """
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_tetraR4R3].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _tetraR3R2(cls, spins):
        """
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_tetraR3R2].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _tetraR4R2(cls, spins):
        """
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_tetraR4R2].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _pentaR3(cls, spins):
        """
        五体クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_pentaR3].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _pentaR4(cls, spins):
        """
        五体クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_pentaR4].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _hex(cls, spins):
        """
        六体クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins.prod(axis=-1).mean(axis=-1)

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
        clusters = [cls._null, cls._point,
                    cls._pairR1, cls._pairR2, cls._pairR3, cls._pairR4,
                    cls._triR2R1, cls._triR3R1, cls._triR3R2, cls._triR4R2,
                    cls._tetraR3R1, cls._tetraR4R3, cls._tetraR3R2,
                    cls._tetraR4R2, cls._pentaR3, cls._pentaR4, cls._hex]
        return np.r_[[c(spins) for c in clusters]].T

    @classmethod
    def print_matrix(cls):
        """
        AAAAAA, AAAAAB, AAAABB, AAABAB, AAABBB,
        ABABAB, AABBBB, ABABBB, ABBBBB, BBBBBBの10種を縦軸に取った
        相関関数xiの行列を表示
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
    FCT構造のTetrahedronの相関関数
    N*8*4のmatrixでスピン情報を受け取る
    4つのスピン配列は[R0, R2, R2, R1]の順
    spin flipによる相関関数の変化量が計算ができるよう各点毎の相関関数を
    return
    """
    _idx_pairR1 = np.array([[0, 3], [1, 2]])
    _idx_pairR2 = np.array([[0, 1], [0, 2], [1, 3], [2, 3]])
    _idx_tri = np.array(
        [[x for x in range(4) if x not in [y]] for y in range(4)])

    @staticmethod
    def _null(spins):
        """
        nullクラスター
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        return np.ones(spins.shape[0:-2])

    @staticmethod
    def _point(spins):
        """
        点クラスターの相関関数
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        return spins.mean(-1).mean(-1)

    @classmethod
    def _pairR1(cls, spins):
        """
        第一隣接の二体クラスターの相関関数
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_pairR1].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _pairR2(cls, spins):
        """
        第一隣接の二体クラスターの相関関数
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_pairR2].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _tri(cls, spins):
        """
        第一隣接正三角形クラスターの相関関数
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_tri].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @staticmethod
    def _tetra(spins):
        """
        四面体クラスターの相関関数
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        return spins.prod(axis=-1).mean(axis=-1)


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
            xi: N×6個のxi (各サイト毎の値を出力)
        """
        clusters = [cls._null, cls._point, cls._pairR1, cls._pairR2,
                    cls._tri, cls._tetra]
        return np.r_[[c(spins) for c in clusters]].T

    @classmethod
    def print_matrix(cls):
        """
        AAAA, AAAB, AABB, ABBA, ABBB, BBBBの6つを縦軸に取った相関関数xiの行列を表示
        """
        AAAA = np.array([1, 1, 1, 1] * 8).reshape(1, 8, 4)
        AAAB = np.array([1, 1, 1, -1] * 8).reshape(1, 8, 4)
        ABBA = np.array([1, -1, -1, 1] * 8).reshape(1, 8, 4)
        AABB = np.array([1, 1, -1, -1] * 8).reshape(1, 8, 4)
        ABBB = np.array([1, -1, -1, -1] * 8).reshape(1, 8, 4)
        BBBB = np.array([-1, -1, -1, -1] * 8).reshape(1, 8, 4)
        spins_array = np.array([AAAA, AAAB, ABBA, AABB, ABBB, BBBB])
        xis = []
        print("Tetrahedron")
        for spins in spins_array:
            xi = ["{0:.3f}".format(x) for x in cls.get_averaged_xi(spins)]
            print(" ".join(xi))
            xis.append(xi)
        return xis

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
        """
        randomに配置
        """
        spins = [int(2*round(random.uniform(conc/2., 0.5+conc/2.))-1)
                 for i in range(4)]
        return cls(spins)

    @classmethod
    def L10(cls, _):
        """
        L10構造
        """
        # spins = [-1, -1, 1, 1]
        spins = [0, 1, 1, 0]
        return cls(spins)

    @classmethod
    def L12_A(cls, _):
        """
        A3B L12構造
        """
        spins = [0, 1, 1, 1]
        return cls(spins)

    @classmethod
    def L12_B(cls, _):
        """
        AB3 L12構造
        """
        # spins = [-1, -1, -1, 1]
        spins = [0, 0, 0, 1]
        return cls(spins)

    def __str__(self):
        return str(self.spins)

    def __getitem__(self, index):
        return self.spins[index]


class FaceCenterTetragonal():
    """
    np.array object of 3dim×4sites
    size: cellの大きさ (size**3)
    arrange: 'random', 'L12', 'L10'を指定
    conc: A, B 二元系におけるAの濃度の初期値 (random以外では不使用)
    other variables
        R0: 原点
        R1, R2, R3, R4: 原点に対する第一,第二,第三,第四隣接原子のサイト位置
                        座標位置の確認はprint_coordinateを使う
        TETRA: 四面体クラスターの組み合わせ
                   原点 + R2から2種 + R1から1種
        OCTA1: 八面体クラスターの組み合わせ 4種
                   [R3×2, R3×2, R4×2]の順番
        OCTA2: 八面体クラスターの組み合わせ 2種
                   [R4×2, R3×2, R3×2]の順番
        TRANS: 座標変換に用いる行列
               原点を別の副格子に移して他の座標を参照する際に使用
               [新しい原点の副格子, 参照する副格子]として記述している
               例: 新しい原点をsite3からsite2を参照する場合
                   + [1, 0, 0, -1]した座標がそのサイトになる
        COORDS: 全ての格子点の座標
        COORDS_ID0: site_id = 0のみの座標

    """
    R0 = [0, 0, 0, 0]
    R1 = [[0, 0, 0, 3], [0, -1, 0, 3], [-1, -1, 0, 3], [-1, 0, 0, 3]]
    R2 = [[0, 0, 0, 1], [0, 0, 0, 2], [0, -1, 0, 1], [-1, 0, 0, 2],
          [0, 0, -1, 1], [0, 0, -1, 2], [0, -1, -1, 1], [-1, 0, -1, 2]]
    R3 = [[1, 0, 0, 0], [-1, 0, 0, 0], [0, 1, 0, 0], [0, -1, 0, 0]]
    R4 = [[0, 0, 1, 0], [0, 0, -1, 0]]

    TETRA = [[R0, R2[0], R2[1], R1[0]], [R0, R2[3], R2[0], R1[3]],
             [R0, R2[2], R2[3], R1[2]], [R0, R2[1], R2[2], R1[1]],
             [R0, R2[4], R2[5], R1[0]], [R0, R2[7], R2[4], R1[3]],
             [R0, R2[6], R2[7], R1[2]], [R0, R2[5], R2[6], R1[1]]]
    OCTA1 = [[R0, R3[0], R1[0], R1[1], R2[1], R2[5]],
             [R0, R3[1], R1[2], R1[3], R2[3], R2[7]],
             [R0, R3[2], R1[0], R1[3], R2[0], R2[4]],
             [R0, R3[3], R1[1], R1[2], R2[2], R2[6]]]
    OCTA2 = [[R0, R4[0], R2[0], R2[2], R2[1], R2[3]],
             [R0, R4[1], R2[4], R2[6], R2[5], R2[7]]]

    TRANS = np.array(
        [[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
         [[0, 0, 0, 1], [0, 1, 1, -1], [0, 0, 1, 1], [0, 1, 0, -1]],
         [[0, 0, 0, 2], [0, 0, 1, 2], [1, 0, 1, -2], [1, 0, 0, -2]],
         [[0, 0, 0, 3], [0, 1, 0, 1], [1, 0, 0, -1], [1, 1, 0, -3]]])

    def __init__(self, size, arrange='random', conc=0.5):
        self.size = size
        mode = {'random': QuadSite.random, 'L10': QuadSite.L10,
                'L12_A': QuadSite.L12_A, 'L12_B': QuadSite.L12_B,}[arrange]
        self.cell = np.array([[[mode(conc).spins for z in range(size)]
                               for y in range(size)] for x in range(size)])
        self.COORDS = np.fromiter(chain.from_iterable(
            product(range(0, self.size), repeat=3)), np.int).reshape(-1, 3)
        self.COORDS_ID0 = np.c_[self.COORDS,
                                np.zeros_like(self.COORDS[:, 0], np.int)]

    def get_octa_xi_glob(self):
        """
        系内全域の相関関数の平均値をreturn
        全てのcoordsにおけるspin配列を作成してget_averaged_xiを利用する
        other variables
            octa1: OCTA1に属する八面体クラスターの格子点
                   (size^3点×6組×4点×4座標)
            octa2: OCTA2に属する八面体クラスターの格子点
                   (size^3点×6組×2点×4座標)
            spins*: spin配列
                    spins2に関しては[R3,R3,R4]の順になるよう一旦、並び替える
        """
        octa1 = (np.array(self.COORDS_ID0).reshape(self.size**3, 1, 1, 4) +
                 self.OCTA1).reshape(self.size**3*6*4, 4)
        octa2 = (np.array(self.COORDS_ID0).reshape(self.size**3, 1, 1, 4) +
                 self.OCTA2).reshape(self.size**3*6*2, 4)
        spins1 = np.r_[[self.from_origin(octa1, i)
                        for i in range(4)]].reshape(self.size**3*4, 24)
        spins2 = np.r_[[self.from_origin(octa2, i)
                        for i in range(4)]].reshape(self.size**3*4, 2, 6)
        spins2 = spins2[..., [2, 3, 4, 5, 0, 1]].reshape(self.size**3*4, 12)
        spinsr = np.c_[spins1, spins2].reshape(self.size**3*4, 6, 6)
        return Octahedron.get_averaged_xi(spinsr)

    def get_tetra_xi_glob(self):
        """
        系内全域の相関関数の平均値をreturn
        全てのcoordsにおけるspin配列を作成してget_averaged_xiを利用する
        other variables
            tetra: 全ての格子点における四面体クラスターのspin配列
                  (size^3点×8組×4点×4座標)
            spins: site id 4種について八面体クラスターspinの6つの組み合わせを
                   求めて行列として結合したもの(size^3*4点×8組×4spin)
        """
        tetra = (np.array(self.COORDS_ID0).reshape(self.size**3, 1, 1, 4) +
                 self.TETRA).reshape(self.size**3*8*4, 4)
        spins = np.r_[[self.from_origin(tetra, i)
                       for i in range(4)]].reshape(self.size**3*4, 8, 4)
        return Tetrahedron.get_averaged_xi(spins)

    def get_TO_xi_glob(self):
        """

        """
        xi_tet = self.get_tetra_xi_glob()
        xi_octa = self.get_octa_xi_glob()
        return np.r_[[xi_octa[0], xi_octa[1], xi_octa[3],
                      xi_octa[2], xi_octa[5], xi_octa[4],
                      xi_octa[6], xi_octa[9], xi_octa[8],
                      xi_octa[7], xi_tet[5], xi_octa[13],
                      xi_octa[12], xi_octa[11], xi_octa[10],
                      xi_octa[15], xi_octa[14], xi_octa[16]]]

    def get_TO_dxi(self, site_id):
        """
        spinがフリップすることによる相関関数の各座標毎の変化量を
        すべての座標に渡って算出する
        Tetrahedron-Octahederon
        return
            size^3 × 11dxi
        other variables
            tetra: 全ての格子点における四面体クラスターのspin配列
                   (size^3点×8組×4点×4座標)
            spins_t: 四面体クラスターspinの8つの組み合わせを
                     求めて行列として結合したもの(size^3点×8組×4spin)
            flipped_t: 先頭のspinがflipしたもの
            octa: 全ての格子点における八面体クラスターのspin配列
                   (size^3点×6組×6点×4座標)
            spins_o: 八面体クラスターspinの6つの組み合わせを
                     求めて行列として結合したもの(size^3点×6組×6spin)
            flipped_o: 先頭のspinがflipしたもの
        """
        tetra = (np.array(self.COORDS_ID0).reshape(self.size**3, 1, 1, 4) +
                 self.TETRA).reshape(self.size**3*8*4, 4)
        spins_t = self.from_origin(tetra, site_id).reshape(self.size**3, 8, 4)
        flipped_t = np.copy(spins_t)
        flipped_t[:, :, 0] *= -1
        dxi_tet = Tetrahedron.get_xi(flipped_t) - Tetrahedron.get_xi(spins_t)

        octa1 = (np.array(self.COORDS_ID0).reshape(self.size**3, 1, 1, 4) +
                 self.OCTA1).reshape(self.size**3*4*6, 4)
        spins_o1 = self.from_origin(octa1, site_id).reshape(self.size**3, 4, 6)
        flipped_o1 = np.copy(spins_o1)
        flipped_o1[:, :, 0] *= -1

        octa2 = (np.array(self.COORDS_ID0).reshape(self.size**3, 1, 1, 4) +
                 self.OCTA2).reshape(self.size**3*2*6, 4)
        spins_o2 = self.from_origin(octa2, site_id).reshape(self.size**3, 2, 6)
        flipped_o2 = np.copy(spins_o2)
        flipped_o2[:, :, 0] *= -1

        spins_o1 = spins_o1.reshape(self.size**3, 24)
        spins_o2 = spins_o2[..., [2, 3, 4, 5, 0, 1]].reshape(self.size**3, 12)
        spins_o = np.c_[spins_o1, spins_o2].reshape(self.size**3, 6, 6)

        flipped_o1 = flipped_o1.reshape(self.size**3, 24)
        flipped_o2 = flipped_o2[...,
                                [2, 3, 4, 5, 0, 1]].reshape(self.size**3, 12)
        flipped_o = np.c_[flipped_o1, flipped_o2].reshape(self.size**3, 6, 6)

        dxi_octa = Octahedron.get_xi(flipped_o) - Octahedron.get_xi(spins_o)
        return np.r_[[dxi_octa[:, 0], dxi_octa[:, 1], dxi_octa[:, 3],
                      dxi_octa[:, 2], dxi_octa[:, 5], dxi_octa[:, 4],
                      dxi_octa[:, 6], dxi_octa[:, 9], dxi_octa[:, 8],
                      dxi_octa[:, 7], dxi_tet[:, 5], dxi_octa[:, 13],
                      dxi_octa[:, 12], dxi_octa[:, 11], dxi_octa[:, 10],
                      dxi_octa[:, 15], dxi_octa[:, 14], dxi_octa[:, 16]]].T

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
        return self.cell[tuple(rpos.T)]

    @staticmethod
    def conv2real_coordinates(coords):
        """
        実空間の座標に変換
        """
        trans = np.array([[0, 0, 0], [0, 0.5, 0.5],
                          [0.5, 0, 0.5], [0.5, 0.5, 0]])
        real_coords = coords[:, :3] + trans[coords[:, 3]]
        return real_coords

    @classmethod
    def print_coordinates(cls):
        """
        R1, R2, TETRA, OCTAのcoordinatesをprint
        確認用のmethod
        """
        print("R1")
        print(cls.conv2real_coordinates(np.array(cls.R1)))
        print("R2")
        print(cls.conv2real_coordinates(np.array(cls.R2)))
        print("TETRA")
        print(cls.conv2real_coordinates(
            np.array(cls.TETRA).reshape(-1, 4)).reshape(-1, 4, 3))
        print("OCTA1")
        print(cls.conv2real_coordinates(
            np.array(cls.OCTA1).reshape(-1, 4)).reshape(-1, 6, 3))
        print("OCTA2")
        print(cls.conv2real_coordinates(
            np.array(cls.OCTA2).reshape(-1, 4)).reshape(-1, 6, 3))

    def get_conc(self, s_direction):
        return (self.cell == s_direction).sum() / self.size ** 3 / 4

    def make_poscar(self, fname):
        """
        self.cellをPOSCARのformatで出力する
        attributes
            fname: 出力のファイル名
        other variables
            poss1, poss2: 元素A, Bの座標
            num1, num2: 元素A, Bのtotalの元素数
            lines: 出力
        """
        poss1 = self.conv2real_coordinates(
            np.argwhere(self.cell == 1)) / self.size
        poss2 = self.conv2real_coordinates(
            np.argwhere(self.cell != 1)) / self.size
        num1 = poss1.shape[0]
        num2 = poss2.shape[0]
        lines = "mc\n"
        lines += "{0}\n".format(self.size*2.5)
        lines += "1.00 0 0\n0 1.00 0\n0 0 1.00\n"
        # lines += "Al Cu\n"
        lines += "Cu\n"
        # lines += "{0} {1}\n".format(num1, num2)
        lines += "{0}\n".format(num2)
        lines += "Direct\n"
        # for site in poss1:
        #     lines += " ".join([str(x) for x in site]) + "\n"
        for site in poss2:
            lines += " ".join([str(x) for x in site]) + "\n"
        with open(fname, 'w') as wfile:
            wfile.write(lines)

    def save_cell(self, title):
        """
        cellをpickleデータに保存する
        """
        fname = title + ".pickle"
        with open(fname, 'wb') as wbfile:
            pickle.dump(self, wbfile)

    @staticmethod
    def load_cell(pickle_path):
        """
        cellをpickleデータからロードする
        """
        with open(pickle_path, 'rb') as rbfile:
            cell = pickle.load(rbfile)
            return cell



if __name__ == '__main__':
    main()


