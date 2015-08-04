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
    cell = FaceCenterCubic(10, arrange='L10')
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
    四面体と異なってサイトの位置関係で相関関数が異なる
    spinsは(0,1) (2,3) (4,5)が第二隣接の組み合わせとして入力する
    spinsを第二隣接同士の各組み合わせで平均をとって[-1/0/1]*3のlistに
    変換して処理した方が早いかもしれない
    その場合、パターンとして3^3=27通りを分類分けすることになる
    other variables:
        _idx_pairR2: R2同士のpairを組にして並べたindex
        _idx_pairR1: R1のpairの組み合わせを示すindex
                     pair_R2を二組選んで、それらの直積をとる
        _idx_triR1: R1のみからなるtriangle clusterのindex
                    pair_R2の直積をとる
        _idx_triR2: R2を含むtriangle clusterのindex
                    pair_R2を二組選びその中から三つを抽出する組み合わせ
    """
    _idx_pairR2 = np.array(range(6)).reshape(3, 2)
    _idx_pairR1 = np.array(
        [list(product(*x)) for x in
         combinations(_idx_pairR2, 2)]).reshape(-1, 2)
    _idx_triR1 = np.fromiter(chain.from_iterable(
        product(*_idx_pairR2)), np.int).reshape(-1, 3)
    _idx_triR2 = np.array(
        [list(combinations(chain.from_iterable(x), 3))
         for x in combinations(_idx_pairR2, 2)]).reshape(-1, 3)
    _idx_tetra1R2 = np.array(
        [[x for x in range(6) if x not in y] for y in _idx_pairR1])
    _idx_tetra2R2 = np.array(
        [[x for x in range(6) if x not in y] for y in _idx_pairR2])
    _idx_penta = np.array(
        [[x for x in range(6) if x not in [y]] for y in range(6)])

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
    def _triR1(cls, spins):
        """
        第一隣接のみの正三角形クラスターの相関関数
        R2のpair以外から選ぶとR1の組み合わせになる
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_triR1].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _tri2R(cls, spins):
        """
        第二隣接を一つ含む二等辺三角形クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_triR2].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _tetra1R2(cls, spins):
        """
        四面体クラスターの相関関数
        spin = (1, -1)の場合はpair * hexで求まるが
        spin = (1, 0)の場合に使えなくなるので直接求める
        3種のR2 pairから1組選ぶ
        残り2組のR2 pairから1つずつ選ぶ
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_tetra1R2].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _tetra2R2(cls, spins):
        """
        正方形クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_tetra2R2].prod(
            axis=-1).mean(axis=-1).mean(axis=-1)

    @classmethod
    def _penta(cls, spins):
        """
        五体クラスターの相関関数
        attributes
            spins: N×6×6
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_penta].prod(
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
        clusters = [cls._null, cls._point, cls._pairR1, cls._pairR2, cls._triR1,
                    cls._tri2R, cls._tetra1R2, cls._tetra2R2, cls._penta, cls._hex]
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
    _idx_pair = np.fromiter(chain.from_iterable(
        combinations(range(4), 2)), np.int).reshape(-1, 2)
    _idx_tri = np.array(
        [[x for x in range(4) if x not in [y]] for y in range(4)])
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
    def _pair(cls, spins):
        """
        第一隣接の二体クラスターの相関関数
        attributes
            spins: N×8×4
        return
            xi: N個のxi
        """
        return spins[..., cls._idx_pair].prod(
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
        spins = [-1, -1, 1, 1]
        return cls(spins)

    @classmethod
    def L12_A(cls, _):
        """
        A3B L12構造
        """
        spins = [-1, 1, 1, 1]
        return cls(spins)

    @classmethod
    def L12_B(cls, _):
        """
        AB3 L12構造
        """
        spins = [-1, -1, -1, 1]
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

    TRANS = [[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
             [[0, 0, 0, 1], [0, 1, 1, -1], [0, 0, 1, 1], [0, 1, 0, -1]],
             [[0, 0, 0, 2], [0, 0, 1, 2], [1, 0, 1, -2], [1, 0, 0, -2]],
             [[0, 0, 0, 3], [0, 1, 0, 1], [1, 0, 0, -1], [1, 1, 0, -3]]]

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


class FaceCenterCubic(object):
    """
    np.array object of 3dim×4sites
    size: cellの大きさ (size**3)
    arrange: 'random', 'L12', 'L10'を指定
    conc: A, B 二元系におけるAの濃度の初期値 (random以外では不使用)
    other variables
        R0: 原点
        R1, R2: 原点に対する第一,第二隣接原子のサイト位置
                座標位置の確認はprint_coordinateを使う
        TETRA: 四面体クラスターの組み合わせ
                   原点 + R2から2種 + R1から1種
        OCTA: 八面体クラスターの組み合わせ 4種
                   [R2×2, R2×2, R2×2]の組み合わせで並んでいる
        TRANS: 座標変換に用いる行列
               原点を別の副格子に移して他の座標を参照する際に使用
               [新しい原点の副格子, 参照する副格子]として記述している
               例: 新しい原点をsite3からsite2を参照する場合
                   + [1, 0, 0, -1]した座標がそのサイトになる
        COORDS: 全ての格子点の座標
        COORDS_ID0: site_id = 0のみの座標
    """
    R0 = [0, 0, 0, 0]
    R1 = [[0, 0, 0, 1], [0, 0, 0, 2], [0, -1, 0, 1], [-1, 0, 0, 2],
          [0, 0, 0, 3], [0, -1, 0, 3], [-1, -1, 0, 3], [-1, 0, 0, 3],
          [0, 0, -1, 1], [0, 0, -1, 2], [0, -1, -1, 1], [-1, 0, -1, 2]]
    R2 = [[1, 0, 0, 0], [-1, 0, 0, 0],
          [0, 1, 0, 0], [0, -1, 0, 0],
          [0, 0, 1, 0], [0, 0, -1, 0]]

    TETRA = [[R0, R1[0], R1[1], R1[4]], [R0, R1[3], R1[0], R1[7]],
             [R0, R1[2], R1[3], R1[6]], [R0, R1[1], R1[2], R1[5]],
             [R0, R1[8], R1[9], R1[4]], [R0, R1[11], R1[8], R1[7]],
             [R0, R1[10], R1[11], R1[6]], [R0, R1[9], R1[10], R1[5]]]
    OCTA = [[R0, R2[0], R1[1], R1[9], R1[4], R1[5]],
            [R0, R2[1], R1[3], R1[11], R1[6], R1[7]],
            [R0, R2[2], R1[0], R1[8], R1[4], R1[7]],
            [R0, R2[3], R1[2], R1[10], R1[5], R1[6]],
            [R0, R2[4], R1[0], R1[2], R1[1], R1[3]],
            [R0, R2[5], R1[8], R1[10], R1[9], R1[11]]]

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

    def make_poscar_order_furface(self, fname, ecis, delta=0.032):
        """
        POSCARを作成する
        規則相の界面をflipのprobabilityから判定
        界面に第三の原子を置いて出力する
        ざっと使ってみたがあまりうまくいかない
        A, B双方からのflipの差が小さくなるよう
        chemical potentialの値を設定する必要がある
        attributes
            fname: 出力のファイル名
            ecis: エネルギー算出の為のeci
            delta: 判定に用いるエネルギーの閾値
        other variables
            delta_xi = フリップに依る相関関数の変化 (TO近似)
            de: フリップに依るエネルギー差
            prob: 遷移確率
            poss: 元素A, B, 界面の座標
            num: 元素A, B, 界面のtotalの元素数
            lines: 出力
        """
        delta_xi = np.hstack([self.get_TO_dxi(x)
                              for x in range(4)]).reshape(-1, 4, 11)
        de = (delta_xi * ecis).sum(axis=2).reshape(
            self.size, self.size, self.size, 4)
        prob = [(self.cell == 1) * (de > delta),
                (self.cell != 1) * (de > delta),
                de <= delta]
        poss = [self.conv2real_coordinates(np.argwhere(x))/self.size
                for x in prob]
        num = [x.shape[0] for x in poss]
        lines = "mc\n"
        lines += "{0}\n".format(self.size*2.5)
        lines += "1.00 0 0\n0 1.00 0\n0 0 1.00\n"
        lines += "Al Cu Ni\n"
        lines += "{0} {1} {2}\n".format(*num)
        lines += "Direct\n"
        for pos in poss:
            for site in pos:
                lines += " ".join([str(x) for x in site]) + "\n"
        with open(fname, 'w') as wfile:
            wfile.write(lines)

    @staticmethod
    def conv2real_coordinates(coords):
        """
        get coordinate in real space
        サイトidに並進対象ベクトルを足してreturn
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
        print("OCTA")
        print(cls.conv2real_coordinates(
            np.array(cls.OCTA).reshape(-1, 4)).reshape(-1, 6, 3))

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

    def get_orderparam(self, cell_range):
        """
        order parameterをreturnする
        order parameterの評価方法は澤田さんのメモを利用
        attributes
            cell_range: 評価するcellの範囲 2*cell_range+1を取る
        return
            (s1, s2, s3): oreder parameters

        local variables
            latt: 2*cell_range+1内の格子点
            idx: 格子点の和が偶数のみのサイトを選び出す
            latt_00: 原点のサイトの組み合わせ
            latt_8: order parameterを計算するための8種のサイト
            p_8: 8種のサイトの占有率
        """
        latt = np.fromiter(chain.from_iterable(
            product(range(2*cell_range+1), repeat=3)), np.int).reshape(-1, 3)
        idx = latt.sum(axis=1) % 2 == 0
        latt_00 = np.c_[latt[idx], [0]*latt[idx].shape[0]]
        trans = [[0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 2], [0, 0, 0, 3],
                 [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 0, 2], [1, 0, 0, 3]]

        latt_8 = [latt_00 + x for x in trans]

        p_8 = [self.cell[tuple(x.T)].mean() for x in latt_8]

        s1 = (p_8[0] + p_8[1] - p_8[2] - p_8[3] +
              p_8[4] + p_8[5] - p_8[6] - p_8[7]) / 8
        s2 = (p_8[0] - p_8[1] + p_8[2] - p_8[3] +
              p_8[4] - p_8[5] - p_8[6] + p_8[7]) / 8
        s3 = (p_8[0] - p_8[1] - p_8[2] + p_8[3] +
              p_8[4] - p_8[5] - p_8[6] + p_8[7]) / 8
        return s1, s2, s3

    def get_TO_xi_glob(self):
        """
        Tetrahedron-Octahedronのグローバルな相関関数を合わせて出力
        nullから八面体までの11種
        get_octa_xi_globとget_tetra_xi_globを呼び出して組み合わせる
        """
        octa = self.get_octa_xi_glob()
        tetra = self.get_tetra_xi_glob()
        return np.array([tetra[0], tetra[1], tetra[2], octa[3], tetra[3],
                         octa[5], tetra[4], octa[6], octa[7], octa[8],
                         octa[9]])

    def get_octa_xi_single(self, coord):
        """
        座標(coord)周りの6つのoctahedronクラスターの相関関数を計算し
        その平均をreturnする
        attributes
            coord: (x, y, z, site_id)
        other variables
            octa: coord 周りの octahedron の座標 (6組×6点×4座標)
            spins: 八面体クラスターspinの6つの組み合わせ (1site×6組×6spin)
        確認用のmethod
        """
        octa = (np.array(coord[:3]+[0]) + self.OCTA).reshape(36, 4)
        spins = self.from_origin(octa, coord[3]).reshape(1, 6, 6)
        return Octahedron.get_averaged_xi(spins)

    def get_octa_xi_glob(self):
        """
        系内全域の相関関数の平均値をreturn
        全てのcoordsにおけるspin配列を作成してget_averaged_xiを利用する
        other variables
            octa: 全ての格子点を含む八面体クラスターのspin配列
                  (size^3点×6組×6点×4座標)
            spins: site id 4種について八面体クラスターspinの6つの組み合わせを
                   求めて行列として結合したもの(size^3*4点×6組×6spin)
        """
        octa = (np.array(self.COORDS_ID0).reshape(self.size**3, 1, 1, 4) +
                self.OCTA).reshape(self.size**3*6*6, 4)
        spins = np.r_[[self.from_origin(octa, i)
                       for i in range(4)]].reshape(self.size**3*4, 6, 6)
        return Octahedron.get_averaged_xi(spins)

    def get_tetra_xi_single(self, coord):
        """
        座標(coord)周りの8つのtetragonalクラスターの相関関数を計算して
        その平均をreturnする
        attributes
            coord: (x, y, z, site_id)
        other variables
            octa: coord 周りの octahedron の座標 (8組×4点×4座標)
            spins: 四面体クラスターspinの8つの組み合わせ (1site×8組×4spin)
        確認用のmethod
        """
        tetra = (np.array(coord[:3]+[0]) + self.TETRA).reshape(32, 4)
        spins = self.from_origin(tetra, coord[3]).reshape(1, 8, 4)
        return Tetrahedron.get_averaged_xi(spins)

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

    def get_tetra_dxi_glob(self):
        """
        spinがフリップすることによる相関関数の変化量をすべての格子点に渡って算出する
        8つのtetrahedronのみを計算
        return
            size^3*4 × 5dxi
        other variables
            tetra: 全ての格子点における四面体クラスターのspin配列
                   (size^3点×8組×4点×4座標)
            spins: site id 4種について八面体クラスターspinの6つの組み合わせを
                   求めて行列として結合したもの(size^3*4点×8組×4spin)
            flipped: 先頭のspinがflipしたもの
        """
        tetra = (np.array(self.COORDS_ID0).reshape(self.size**3, 1, 1, 4) +
                 self.TETRA).reshape(self.size**3*8*4, 4)
        spins = np.r_[[self.from_origin(tetra, i)
                       for i in range(4)]].reshape(self.size**3*4, 8, 4)
        flipped = np.copy(spins)
        flipped[:, :, 0] = spins[:, :, 0] * -1
        return (Tetrahedron.get_xi(flipped) -
                Tetrahedron.get_xi(spins))

    def get_tetra_dxi(self, site_id):
        """
        spinがフリップすることによる相関関数の各site_id毎の変化量を
        すべての座標に渡って算出する
        8つのtetrahedronのみを計算
        return
            size^3 × 5dxi
        other variables
            tetra: 全ての格子点における四面体クラスターのspin配列
                   (size^3点×8組×4点×4座標)
            spins: 八面体クラスターspinの6つの組み合わせを
                   求めて行列として結合したもの(size^3点×8組×4spin)
            flipped: 先頭のspinがflipしたもの
        """
        tetra = (np.array(self.COORDS_ID0).reshape(self.size**3, 1, 1, 4) +
                 self.TETRA).reshape(self.size**3*8*4, 4)
        spins = self.from_origin(tetra, site_id).reshape(self.size**3, 8, 4)
        flipped = np.copy(spins)
        flipped[:, :, 0] *= -1
        return Tetrahedron.get_xi(flipped) - Tetrahedron.get_xi(spins)

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

        octa = (np.array(self.COORDS_ID0).reshape(self.size**3, 1, 1, 4) +
                self.OCTA).reshape(self.size**3*6*6, 4)
        spins_o = self.from_origin(octa, site_id).reshape(self.size**3, 6, 6)
        flipped_o = np.copy(spins_o)
        flipped_o[:, :, 0] *= -1
        dxi_octa = Octahedron.get_xi(flipped_o) - Octahedron.get_xi(spins_o)
        return np.r_[[dxi_tet[:, 0], dxi_tet[:, 1], dxi_tet[:, 2],
                      dxi_octa[:, 3], dxi_tet[:, 3], dxi_octa[:, 5],
                      dxi_tet[:, 4], dxi_octa[:, 6], dxi_octa[:, 7],
                      dxi_octa[:, 8], dxi_octa[:, 9]]].T

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

    def get_conc(self, s_direction):
        """
        スピン(s_direction)の濃度を算出
        """
        return (self.cell == s_direction).sum() / self.size ** 3 / 4

    def __getitem__(self, *pos):
        return self.cell.__getitem__(*pos)

if __name__ == '__main__':
    main()


