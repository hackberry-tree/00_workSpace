#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
from __future__ import division
import math
import pickle
import random
import numpy as np
from itertools import product, chain, combinations, permutations


class MonteCarlo(object):
    """
    ECIs, cell, 温度 T をセットに読み込んで, Monte-Carlo の loop 計算を行う
    kb はボルツマン定数
    default は meV/K
    """
    def __init__(self, cell, T, kb=8.6171e-2):
        self.cell = cell
        self.kb = kb
        if T == 0:
            self.trans_prob = self.trans_prob_test
            self.T = 0
        else:
            self.T = T
            self.trans_prob = self.trans_prob_de
        # self.beta = self.kb * self.T

    @property
    def beta(self):
        return self.kb * self.T


    @staticmethod
    def trans_prob_test(de):
        """
        de が negative の場合必ず flip する
        test 用
        """
        de[de >= 0] = 0
        de[de < 0] = 0.01
        return de

    def trans_prob_de(self, de):
        """
        エネルギー差の行列から遷移確率の行列を返す
        """
        w = np.exp(-de/self.beta)
        return w / (1+w)

    def loop_fcc_grand(self, steps):
        """
        fcc cell の spin-flip loop
        grand canonical ensemble
        """
        # energy = [[self.cell.get_conc(), self.cell.get_energy()]]
        print("initial concentration")
        print(self.cell.get_conc())
        print("initial energy")
        print(self.cell.get_energy())
        print("start")
        for _ in range(steps):
            de = self.cell.get_flip_de()
            p = self.trans_prob(de)
            rand = np.random.rand(p.size).reshape(
                self.cell.size, self.cell.size, self.cell.size, 4)
            judge = (rand < p)
            self.cell.cell[judge] = 1 - self.cell.cell[judge]
            print(judge.sum())
            print(self.cell.get_energy())
            print(self.cell.get_conc())
            # energy.append([self.cell.get_conc(), self.cell.get_energy()])
        # return energy

    def loop_fcc_micro(self, steps):
        """
        fcc cell の spin flip loop
        micro canonical ensemble
        """
        print("initial concentration")
        print(self.cell.get_conc())
        print("initial energy")
        print(self.cell.get_energy())
        print("start")
        # energy = [[0, self.cell.get_energy()]]

        for _ in range(steps):
            de = self.cell.get_flip_de()
            m2p = self.cell.cell == 0
            p2m = self.cell.cell == 1
            is_m2p_major = m2p.sum() > p2m.sum()
            major = {True: m2p, False: p2m}[is_m2p_major]
            minor = {True: p2m, False: m2p}[is_m2p_major]
            pair = np.random.choice(
                range(major.sum()), minor.sum(), replace=False)
            delta = de[tuple(np.array(np.where(major))[:, pair])] + \
                de[tuple(np.array(np.where(minor)))]
            p = self.trans_prob(delta)
            rand = np.random.rand(p.size)
            judge = (rand < p)
            jmat = major * False
            jmat[tuple(np.array(np.where(major))[:, pair])] = judge
            jmat[tuple(np.array(np.where(minor)))] = judge
            self.cell.cell[jmat] = 1 - self.cell.cell[jmat]
            print(self.cell.get_energy())

            a = judge.sum()
            # energy.append([a, self.cell.get_energy()])
            print(a)
        # return energy

    def __select_pair_sites(self):
        """
        s = 0/1 の site の pair をランダムに選択する関数を生成する
        """
        num_s0 = range((self.cell.cell == 0).sum())
        num_s1 = range((self.cell.cell == 1).sum())
        def select_pair_sites():
            """
            num_s0, s1 を memo 化しておく
            """
            s0 = self.cell.cell == 0
            site_s0 = (
                np.array(np.where(s0)))[:, np.random.choice(
                    num_s0, 1, replace=False)[0]]
            s1 = self.cell.cell == 1
            site_s1 = (
                np.array(np.where(s1)))[:, np.random.choice(
                    num_s1, 1, replace=False)[0]]
            return site_s0, site_s1
        return select_pair_sites

    def __loop_fcc_micro_single(self, steps, func):
        """
        fcc cell の spin flip loop
        micro canonical ensemble
        1 点づつ判定する
        エネルギー差を計算する func を指定する

        # 01 s=0,1 のサイトをランダムに選択
        # 02 spin を交換した時のエネルギーを計算 -> flip 判定
        # 03 spin の交換処理

        """
        # print("initial concentration")
        # print(self.cell.get_conc())
        # print("initial energy")
        # total_e = self.cell.get_energy()
        # print(total_e)
        # print("start")

        # 01
        select_pair_sites = self.__select_pair_sites()
        # size = self.cell.size**3*4
        delta_e = 0
        size = self.cell.size ** 3 * 4
        # initial = self.cell.get_energy()
        for _ in range(steps):
            # 01
            site_s0, site_s1 = select_pair_sites()
            # 02
            de = func(site_s0, site_s1)
            p = float(self.trans_prob(de))
            # print(p)
            rand = float(np.random.rand(1))
            s = {False: [0, 1], True: [1, 0]}[rand < p]

            # 03
            self.cell.cell[tuple(site_s0.T)] = s[0]
            self.cell.cell[tuple(site_s1.T)] = s[1]
            out_de = {False: 0, True: de}[rand < p]
            delta_e += out_de
        return delta_e
        # print(delta_e/size)
        # final = self.cell.get_energy()
        # print(final - initial)

    def loop_fcc_micro_single(self, steps):
        """
        normal size の cell 用
        """
        return self.__loop_fcc_micro_single(steps, self.cell.get_exchange_de)

    def loop_fcc_micro_single_small(self, steps):
        """
        NaCl 型 cell の spin flip loop
        micro canonical ensemble
        1 点づつ判定
        小さな cell 用の計算 (cell のサイズが相関関数のサイズより小さい場合)
        詳しくは get_exchange_de_small を参照
        """
        self.__loop_fcc_micro_single(steps, self.cell.get_exchange_de_small)

    def loop_fcc_micro_prio_de(self, steps):
        """
        NaCl 型 cell の spin flip loop
        micro canonical ensemble
        エネルギーが低いものを優先して flip
        書きかけ
        """
        print("initial concentration")
        print(self.cell.get_conc())
        print("initial energy")
        print(self.cell.get_energy())
        print("start")
        # energy = [[0, self.cell.get_energy()]]

        for _ in range(steps):
            de = self.cell.get_flip_de()
            m2p = self.cell.cell == 0
            p2m = self.cell.cell == 1
            # is_m2p_major = m2p.sum() > p2m.sum()
            # major = {True: m2p, False: p2m}[is_m2p_major]
            # minor = {True: p2m, False: m2p}[is_m2p_major]
            # print(np.sort(de[m2p]))
            # print((de == np.sort(de[p2m])[0]).sum())
            rand = np.random.choice(range(5), 1)
            rand2 = np.random.choice(range(5), 1)
            rand = 0
            rand2 = 0
            self.cell.cell[de == np.sort(de[p2m])[rand2]] = 1
            self.cell.cell[de == np.sort(de[p2m])[rand]] = 0
            print(self.cell.get_energy())

    def loop_nacl_grand(self, steps):
        """
        NaCl 型 cell の spin flip loop
        grand canonical ensemble
        """
        # print("initial concentration")
        # print(self.cell.get_conc())
        # print("initial energy")
        # print(self.cell.get_energy())
        # print("start")

        for _ in range(steps):
            de = self.cell.get_flip_de_sub()
            p = self.trans_prob(de)
            rand = np.random.rand(p.size)
            judge = (rand < p).reshape(4, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 4)
            self.cell.cell[..., 0:4][judge] = 1 - \
                self.cell.cell[..., 0:4][judge]

            de = self.cell.get_flip_de_int()
            p = self.trans_prob(de)
            rand = np.random.rand(p.size)
            judge = (rand < p).reshape(4, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 4)
            self.cell.cell[..., 4:8][judge] = 1 - \
                self.cell.cell[..., 4:8][judge]
        #     print(self.cell.get_energy())
        # print(self.cell.get_conc())

    def loop_flip_nacl_fix_conc(self, steps):
        """
        NaCl 型 cell の spin flip loop
        spin flip の数を等しくする
        このやり方は好ましくない
        ToDo: log を残して消す
        """
        print("initial concentration")
        print(self.cell.get_conc())
        print("initial energy")
        print(self.cell.get_energy())
        print("start")

        for _ in range(steps):
            de = self.cell.get_flip_de_sub()
            p = self.trans_prob(de)
            rand = np.random.rand(p.size)
            judge = (rand < p).reshape(4, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 4)
            m2p = judge * (self.cell.cell[..., 0:4] == 0)
            p2m = judge * (self.cell.cell[..., 0:4] == 1)
            is_m2p_major = m2p.sum() > p2m.sum()
            major = {True: m2p, False: p2m}[is_m2p_major]
            minor = {True: p2m, False: m2p}[is_m2p_major]
            diff = major.sum() - minor.sum()
            if diff != 0:
                fix = np.random.choice(
                    range(major.sum()), diff, replace=False)
                major[tuple(np.array(np.where(major))[:, fix])] = False
            t = major + minor
            print(t.sum(), diff)
            self.cell.cell[..., 0:4][t] = 1 \
                - self.cell.cell[..., 0:4][t]

            de = self.cell.get_flip_de_int()
            p = self.trans_prob(de)
            rand = np.random.rand(p.size)
            judge = (rand < p).reshape(4, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 4)
            m2p = judge * (self.cell.cell[..., 4:8] == 0)
            p2m = judge * (self.cell.cell[..., 4:8] == 1)
            is_m2p_major = m2p.sum() > p2m.sum()
            major = {True: m2p, False: p2m}[is_m2p_major]
            minor = {True: p2m, False: m2p}[is_m2p_major]
            diff = major.sum() - minor.sum()
            if diff != 0:
                fix = np.random.choice(
                    range(major.sum()), diff, replace=False)
                major[tuple(np.array(np.where(major))[:, fix])] = False
            t = major + minor
            print(t.sum(), diff)
            self.cell.cell[..., 4:8][t] = 1 - \
                self.cell.cell[..., 4:8][t]

            print(self.cell.get_energy())
        print(self.cell.get_conc())

    def loop_nacl_micro(self, steps):
        """
        NaCl 型 cell の spin flip loop
        micro canonical ensemble
        """
        print("initial concentration")
        print(self.cell.get_conc())
        print("initial energy")
        print(self.cell.get_energy())
        print("start")
        energy = [[0, 0, self.cell.get_energy()]]

        for _ in range(steps):
            de = self.cell.get_flip_de_sub().reshape(4, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 4)
            m2p = self.cell.cell[..., 0:4] == 0
            p2m = self.cell.cell[..., 0:4] == 1
            is_m2p_major = m2p.sum() > p2m.sum()
            major = {True: m2p, False: p2m}[is_m2p_major]
            minor = {True: p2m, False: m2p}[is_m2p_major]
            pair = np.random.choice(
                range(major.sum()), minor.sum(), replace=False)
            delta = de[tuple(np.array(np.where(major))[:, pair])] + \
                de[tuple(np.array(np.where(minor)))]
            p = self.trans_prob(delta)
            rand = np.random.rand(p.size)
            judge = (rand < p)
            jmat = major * False
            jmat[tuple(np.array(np.where(major))[:, pair])] = judge
            jmat[tuple(np.array(np.where(minor)))] = judge
            self.cell.cell[..., 0:4][jmat] = 1 - self.cell.cell[..., 0:4][jmat]

            a = judge.sum()

            de = self.cell.get_flip_de_int().reshape(4, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 4)
            m2p = self.cell.cell[..., 4:8] == 0
            p2m = self.cell.cell[..., 4:8] == 1
            is_m2p_major = m2p.sum() > p2m.sum()
            major = {True: m2p, False: p2m}[is_m2p_major]
            minor = {True: p2m, False: m2p}[is_m2p_major]
            pair = np.random.choice(
                range(major.sum()), minor.sum(), replace=False)
            delta = de[tuple(np.array(np.where(major))[:, pair])] + \
                de[tuple(np.array(np.where(minor)))]
            p = self.trans_prob(delta)
            rand = np.random.rand(p.size)
            judge = (rand < p)
            jmat = major * False
            jmat[tuple(np.array(np.where(major))[:, pair])] = judge
            jmat[tuple(np.array(np.where(minor)))] = judge
            self.cell.cell[...,4:8][jmat] = 1 - self.cell.cell[..., 4:8][jmat]

            b = judge.sum()
            energy.append([a, b, self.cell.get_energy()])
            print(self.cell.get_energy())
        return energy


class NaClSite(object):
    """
    NaCl型の 8 点のサイトの情報を記録する
    s-site 4 点 [0:4]
    i-site 4 点 [4:8]
    args:
        spins: spin の配列を記した 8 つの list
    cls variables:
        SPIN: spin の取り方を指定 default は [1, 0]
        SITES: 各サイトの座標
        CONC_INT: i-site の 1-spin の濃度
        CONC_SUB: s-site の 1-spin の濃度
    """
    SPIN = [1, 0]
    _site_s0 = [0.0, 0.0, 0.0]
    _site_s1 = [0.0, 1/2, 1/2]
    _site_s2 = [1/2, 0.0, 1/2]
    _site_s3 = [1/2, 1/2, 0.0]
    _site_i0 = [1/2, 0.0, 0.0]
    _site_i1 = [0.0, 1/2, 0.0]
    _site_i2 = [0.0, 0.0, 1/2]
    _site_i3 = [1/2, 1/2, 1/2]
    SITES = [_site_s0, _site_s1, _site_s2, _site_s3,
             _site_i0, _site_i1, _site_i2, _site_i3]
    CONC_INT = 0.5
    CONC_SUB = 0.5

    def __init__(self, spins):
        self.spins = spins

    @classmethod
    def conv_site2idex(cls, site):
        """
        座標の直接表記を class 内の index 表記に変換して return する
        """
        index = [math.floor(x) for x in site]
        decim = [x - math.floor(x) for x in site]
        index.append({tuple(cls.SITES[i]): i
                      for i in range(len(cls.SITES))}[tuple(decim)])
        return index

    @classmethod
    def random(cls):
        """
        spin を site に random に配置する
        cls.CONC_INT と cls.CONC_SUB で i-, s-site の濃度を指定する
        ToDo: cls.SPIN に対応させる
        """
        i = cls.CONC_INT
        s = cls.CONC_SUB
        spins = [cls.SPIN[int(round(random.uniform(0.5-s/2., 1-s/2.)))]
                 for _ in range(4)] + \
                [cls.SPIN[int(round(random.uniform(0.5-i/2., 1-i/2.)))]
                 for _ in range(4)]
        return cls(spins)

    @classmethod
    def FCC_u(cls):
        """
        FCC 構造
        s-site が up-spin
        """
        i, j = cls.SPIN
        spins = [i, i, i, i, j, j, j, j]
        return cls(spins)

    @classmethod
    def FCC_d(cls):
        """
        FCC 構造
        s-site が down-spin
        """
        _, j = cls.SPIN
        spins = [j, j, j, j, j, j, j, j]
        return cls(spins)

    @classmethod
    def L10(cls):
        """
        L10 構造
        """
        i, j = cls.SPIN
        spins = [i, i, j, j, j, j, j, j]
        return cls(spins)

    @classmethod
    def L12_A(cls):
        """
        A3B L12 構造
        """
        i, j = cls.SPIN
        spins = [i, i, i, j, j, j, j, j]
        return cls(spins)

    @classmethod
    def NaCl_uu(cls):
        """
        NaCl 型構造
        i-, s-site ともに up-spin
        """
        i, _ = cls.SPIN
        spins = [i, i, i, i, i, i, i, i]
        return cls(spins)

    @classmethod
    def NaCl_du(cls):
        """
        NaCl型構造
        i-site down-spin, s-site up-spin
        """
        i, j = cls.SPIN
        spins = [j, j, j, j, i, i, i, i]
        return cls(spins)


    def __str__(self):
        return str(self.spins)

    def __getitem__(self, index):
        return self.spins[index]


class NaClXtal(object):
    """
    NaCl 型の結晶格子を作成する
    各サイトは index 表記になっているため、座標の参照の仕方が面倒
    4 番目のコラムが unit cell 内の site の index

    0 以外の site から参照する際に self.TRANS を利用する
    例として、site-1 から site-2 を参照させる場合
    self.TRANS[1, 2] を index に足し合わせる
    (直感的に分かりづらいのが問題)

    args:
        size: cell のsize size^3*8sites を作成する
        arrange: 原子配列の初期値
                 選択肢は mode を参照
        conc_sub, conc_int: up-spinの組成比の初期値
                            arrange=random のみ使用
        ecis: ECIs
              format は dict形式 (get_clus.pyから作成する)
        ToDo: get_clus.py をこの module に組み込む
    variables:
        nodup_int_ecis: i-site の ecis で s-site と重複していないもののlist
                        total energy を算出する際に使用
                        これをやらないと重複して相関関数を数えてしまうので、
                        total energy が変わってしまう
    """
    def __init__(self, size, arrange='random',
                 conc_sub=0.5, conc_int=0.5, ecis=None):
        self.size = size
        NaClSite.CONC_INT = conc_int
        NaClSite.CONC_SUB = conc_sub
        mode = {'random': NaClSite.random, 'L10': NaClSite.L10,
                'FCC_u': NaClSite.FCC_u, 'NaCl_uu': NaClSite.NaCl_uu,
                'FCC_d': NaClSite.FCC_d, 'NaCl_du': NaClSite.NaCl_du,}[arrange]
        self.cell = np.array([[[mode().spins for z in range(size)]
                               for y in range(size)] for x in range(size)])

        self.COORDS = np.fromiter(chain.from_iterable(
            product(range(0, self.size), repeat=3)), np.int).reshape(-1, 3)
        self.COORDS_ID0 = np.c_[self.COORDS, np.zeros_like(
            self.COORDS[:, 0], np.int)].reshape(-1, 1, 1, 4)
        self.TRANS = self.get_trans()

        if ecis:
            self.int_clusters = ecis['int_clusters']
            self.sub_clusters = ecis['sub_clusters']
            self.int_ecis = ecis['int_ecis']
            self.sub_ecis = ecis['sub_ecis']
            self.nodup_int_ecis = []
            for e in self.int_ecis:
                if e in self.sub_ecis:
                    self.nodup_int_ecis.append(0)
                else:
                    self.nodup_int_ecis.append(e)
            self.int_idxs = self.get_clus_idxs_int(self.int_clusters)
            self.sub_idxs = self.get_clus_idxs_sub(self.sub_clusters)

    def save_cell(self, title):
        """
        cellをpickleデータに保存する
        """
        fname = title + ".pickle"
        with open(fname, 'wb') as wbfile:
            pickle.dump(self, wbfile)

    @staticmethod
    def get_trans():
        """
        相対座標へ変換するための matrix を return する
        """
        sites = np.array(NaClSite.SITES)
        size = sites.shape[0]
        tmp = sites.reshape(size, 1, 3) + sites.reshape(1, size, 3)
        index = np.array([[NaClSite.conv_site2idex(x) for x in y] for y in tmp])
        obj = np.c_[
            np.zeros_like(np.arange(size*3).reshape(size, 3)), np.arange(size)]
        return index - obj

    @classmethod
    def from_pickle_ecis(cls, fname, size=10, arrange='random',
                         conc_sub=0.5, conc_int=0.5):
        """
        ecis や cluster の情報が入った pickle data から cls を生成
        load する pickle data の format
        dict{sub_clusters: [n * cluster * symm * 4 座標]
             sub_ecis: [n * ecis]
             int_clusters: ...
             int_ecis: ...}
        """
        with open(fname, 'rb') as rbfile:
            ecis_data = pickle.load(rbfile)
        return cls(size, arrange, conc_sub, conc_int, ecis_data)

    def get_xi_sub(self):
        """
        置換型サイトの相関関数を出力
        """
        for idxs in self.sub_idxs:
            spins = self.cell[idxs].T
            xi = spins.prod(axis=-1).mean(axis=-1)
            print(xi.mean())

    def get_xi_int(self):
        """
        侵入型サイトの相関関数を出力
        """
        for idxs in self.int_idxs:
            spins = self.cell[idxs].T
            xi = spins.prod(axis=-1).mean(axis=-1)
            print(xi.mean())


    def get_energy(self):
        """
        系全体のエネルギーを単位元素あたりに換算して return
        """
        ene = np.array(0)
        for idxs, eci in zip(self.sub_idxs, self.sub_ecis):
            spins = self.cell[idxs].T
            xi = spins.prod(axis=-1).mean(axis=-1)
            ene = xi * eci + ene
        for idxs, eci in zip(self.int_idxs, self.nodup_int_ecis):
            spins = self.cell[idxs].T
            xi = spins.prod(axis=-1).mean(axis=-1)
            ene = xi * eci + ene
        return ene.mean()

    def get_flip_de_sub(self):
        """
        置換型サイトの spin flip によるエネルギー差をサイト毎に return
        """
        de = 0
        for idxs, eci in zip(self.sub_idxs, self.sub_ecis):
            spins = self.cell[idxs].T
            xi = spins.prod(axis=-1).mean(axis=-1)
            spins_inv = spins
            spins_inv[..., 0] = 1 - spins[..., 0]
            xi_inv = spins_inv.prod(axis=-1).mean(axis=-1)
            dxi = xi_inv - xi
            de = dxi * eci + de
        return de

    def get_flip_de_int(self):
        """
        侵入型サイトの spin flip によるエネルギー差をサイト毎に return
        """
        de = 0
        for idxs, eci in zip(self.int_idxs, self.int_ecis):
            spins = self.cell[idxs].T
            xi = spins.prod(axis=-1).mean(axis=-1)
            spins_inv = spins
            spins_inv[..., 0] = 1 - spins[..., 0]
            xi_inv = spins_inv.prod(axis=-1).mean(axis=-1)
            dxi = xi_inv - xi
            de = dxi * eci + de
        return de

    def get_clus_idxs_sub(self, clusters):
        """
        fancy index に用いる cluster の index を list 化して return する
        substitutional 用
        """
        indexes = []
        for clus in clusters:
            size = clus.shape
            tmp_pos = []
            for site_id in range(4):
                pos = clus + self.TRANS[site_id][clus[..., 3]] + \
                    self.COORDS_ID0
                tmp_pos.append(pos)
            index = np.array(tmp_pos).reshape(
                4*self.size**3, size[0], size[1], 4)
            index[..., 0:3] %= self.size
            indexes.append(tuple(index.T))
        return indexes

    def get_clus_idxs_int(self, clusters):
        """
        fancy index に用いる cluster の index を list 化して return する
        interstitial 用
        """
        indexes = []
        for clus in clusters:
            size = clus.shape
            tmp_pos = []
            for site_id in range(4, 8):
                pos = clus + self.TRANS[site_id][clus[..., 3]] + \
                    self.COORDS_ID0
                tmp_pos.append(pos)
            index = np.array(tmp_pos).reshape(
                4*self.size**3, size[0], size[1], 4)
            index[..., 0:3] %= self.size
            indexes.append(tuple(index.T))
        return indexes

    def get_conc(self):
        """
        組成比を return
        return:
            {'inter': , 'subs': }
        """
        return {'subs': self.cell[..., 0:4].mean(),
                'inter': self.cell[..., 4:8].mean()}

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
            np.argwhere(self.cell[..., 0:4] == 0), 's') / self.size
        poss2 = self.conv2real_coordinates(
            np.argwhere(self.cell[..., 4:8] == 1), 'i') / self.size
        num1 = poss1.shape[0]
        num2 = poss2.shape[0]
        lines = "mc\n"
        lines += "{0}\n".format(self.size*2.5)
        lines += "1.00 0 0\n0 1.00 0\n0 0 1.00\n"
        lines += "Cr C\n"
        lines += "{0} {1}\n".format(num1, num2)
        lines += "Direct\n"
        for site in poss1:
            lines += " ".join([str(x) for x in site]) + "\n"
        for site in poss2:
            lines += " ".join([str(x) for x in site]) + "\n"
        with open(fname, 'w') as wfile:
            wfile.write(lines)

    @staticmethod
    def conv2real_coordinates(coords, s_or_i):
        """
        実空間の座標に変換
        """
        n = {'s': 0, 'i': 4}[s_or_i]
        trans = np.array(NaClSite.SITES)
        real_coords = coords[:, 0:3] + trans[coords[:, 3] + n]
        return real_coords


class QuadSite(object):
    """
    4点のサイトの情報を記述する
    spins
    sites(不使用、memo)
    """
    SPIN = [1, 0]
    _site_s0 = [0.0, 0.0, 0.0]
    _site_s1 = [0.0, 1/2, 1/2]
    _site_s2 = [1/2, 0.0, 1/2]
    _site_s3 = [1/2, 1/2, 0.0]
    SITES = [_site_s0, _site_s1, _site_s2, _site_s3]
    def __init__(self, spins):
        self.spins = spins

    @classmethod
    def conv_site2idex(cls, site):
        """
        座標の直接表記を class 内の index 表記に変換する
        """
        index = [math.floor(x) for x in site]
        decim = [x - math.floor(x) for x in site]
        index.append({tuple(cls.SITES[i]): i
                      for i in range(len(cls.SITES))}[tuple(decim)])
        return index

    @classmethod
    def random(cls, conc):
        """
        randomに配置
        """
        spins = [int(2*round(random.uniform(conc/2., 0.5+conc/2.))/2)
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

    @classmethod
    def FCC_A(cls, _):
        """
        AB3 L12構造
        """
        # spins = [-1, -1, -1, 1]
        spins = [1, 1, 1, 1]
        return cls(spins)

    def __str__(self):
        return str(self.spins)

    def __getitem__(self, index):
        return self.spins[index]


class FCCXtal(object):
    """
    np.array object of 3dim×4sites
    size: cellの大きさ (size**3)
    arrange: 'random', 'L12', 'L10'を指定
    conc: A, B 二元系におけるAの濃度の初期値 (random以外では不使用)
    other variables
        TRANS: 相対座標
               原点を別の副格子にとった際、他の座標を参照する際に使用
               [新しい原点の副格子, 参照する副格子]として記述している
               例: site3からsite2を参照する場合
                   + [1, 0, 0, -1]した座標がそのサイトになる
        COORDS: 全ての格子点の座標 (site_id は 0 のみ)

    相関関数の計算方法は利便性のために原点[0, 0, 0, 0]には
    必ず spin = 1 を置いて計算する (xi_o とする)
    各サイトの xi_o を予め計算しておいて、
    実際のサイト上の xi_o は サイト上の spin * xi_o で計算する
    spin flip のエネルギー算出の際に計算を簡略化できる

    xi_o 算出のために point cluster は initialize の時点で分離
    (get_clus.py の時点で分離しておく)
    また 各クラスターの原点の座標 [0, 0, 0, 0] も取り除く
    """

    def __init__(self, size, arrange='random', conc=0.5, ecis=None):
        self.size = size
        mode = {'random': QuadSite.random, 'L10': QuadSite.L10,
                'L12_A': QuadSite.L12_A, 'L12_B': QuadSite.L12_B,
                'FCC_A': QuadSite.FCC_A}[arrange]
        self.cell = np.array([[[mode(conc).spins for z in range(size)]
                               for y in range(size)] for x in range(size)])
        self.COORDS = np.fromiter(chain.from_iterable(
            product(range(0, self.size), repeat=3)), np.int).reshape(-1, 3)
        self.COORDS = np.c_[self.COORDS,
                            np.zeros_like(
                                self.COORDS[:, 0], np.int)].reshape(-1, 1, 1, 4)
        self.TRANS = self.get_trans()
        if ecis:
            self.clusters = ecis['clusters']
            self.ecis = ecis['ecis']
            self.eci_point = ecis['eci_point']
            self.idxs = self.get_clus_idxs(self.clusters)

    @classmethod
    def from_pickle_ecis(cls, fname, size=10, arrange='random', conc=0.5):
        """
        ecis や cluster の情報が入った pickle data から cls を生成
        load する pickle data の format は
        dict{clusters: [(n - 2) × degen × (vertex - 1) × 4 座標]
             ecis: [(n - 2) * ecis]
             eci_point: point cluster のみの ecis}
        n-2 は null と point を除いてあることに
        vertex-1 は 原点 [0, 0, 0, 0] を除いてあることに対応する
        """
        with open(fname, 'rb') as rbfile:
            ecis_data = pickle.load(rbfile)
        return cls(size, arrange, conc, ecis_data)

    def get_clus_idxs(self, clusters):
        """
        fancy index に用いる cluster の index を list 化して return する
        index の shape は
        cluster種 × 参照する座標 × クラスターの頂点 ×
            対称操作で複製されたクラスター × すべての座標
        """
        indexes = []
        for clus in clusters:
            size = clus.shape
            tmp_pos = []
            for site_id in range(4):
                pos = clus + self.TRANS[site_id][clus[..., 3]] + \
                    self.COORDS
                tmp_pos.append(pos)
            index = np.array(tmp_pos).reshape(
                4*self.size**3, size[0], size[1], 4)
            index[..., 0:3] %= self.size
            indexes.append(tuple(index.T))
        return indexes

    @staticmethod
    def get_trans():
        """
        相対座標へ変換するための matrix を return する
        """
        sites = np.array(QuadSite.SITES)
        size = sites.shape[0]
        tmp = sites.reshape(size, 1, 3) + sites.reshape(1, size, 3)
        index = np.array([[QuadSite.conv_site2idex(x) for x in y] for y in tmp])
        obj = np.c_[
            np.zeros_like(np.arange(size*3).reshape(size, 3)), np.arange(size)]
        return index - obj

    def _get_site_energy(self):
        """
        siteに仮想的に spin=1 を入れて、その周囲の相関関数を計算
        ecis との積から サイト毎のエネルギーを計算する
        """
        energy = self.eci_point
        for idxs, eci in zip(self.idxs, self.ecis):
            spins = self.cell[idxs].T
            xi = spins.prod(axis=-1).mean(axis=-1)
            energy = (xi * eci) + energy
        return energy.reshape(4, -1).T.reshape(
            self.size, self.size, self.size, 4)

    def get_energy(self):
        """
        total energy を return する
        null の ecis は含まないので注意
        """
        return (self._get_site_energy() * self.cell).mean()

    def _get_site_de(self):
        """
        siteに仮想的に spin=1 を入れて、その周囲の相関関数を計算
        頂点の数をとの積をかけてある
        delta との積から サイト毎の spin flip のエネルギーを計算する
        """
        energy = self.eci_point
        for idxs, eci in zip(self.idxs, self.ecis):
            spins = self.cell[idxs].T
            vertex = spins.shape[2] + 1
            xi = spins.prod(axis=-1).mean(axis=-1)
            energy = (xi * eci) * vertex + energy
        return energy.reshape(4, -1).T.reshape(
            self.size, self.size, self.size, 4)

    def get_flip_de(self):
        """
        spin flip のエネルギー差の matrix を return
        delta は 0 -> 1 の場合 1
                 1 -> 0 の場合 -1 をとる係数
        """
        delta = -2 * self.cell + 1
        return delta * self._get_site_de()

    def get_exchange_de(self, site_s0, site_s1):
        """
        site_s0 は s=0 のサイト s1 は s=1 のサイトを指定
        site_s0, s1 の spin を交換した場合のエネルギー差を return する
        site_s0, s1 の表記は x, y, z, site_id 4 要素の array
        """
        tmp_cell = self.cell.copy()
        tmp_cell[tuple(site_s0.T)] = 1
        tmp_cell[tuple(site_s1.T)] = 0

        energy = 0
        for clus, eci in zip(self.clusters, self.ecis):
            index1 = (clus + self.TRANS[site_s0[3]][clus[..., 3]] +
                      np.r_[site_s0[0:3], [0]])
            index1[..., 0:3] %= self.size
            index2 = (clus + self.TRANS[site_s1[3]][clus[..., 3]] +
                      np.r_[site_s1[0:3], [0]])
            index2[..., 0:3] %= self.size

            # s=1 のサイトの相関関数
            s11 = self.cell[tuple(index2.T)].T
            # s=0 だったところに 1 が入った時の相関関数
            # この時、 s=1 だったところは 0になっているので
            # temp_cell を使う
            s01 = tmp_cell[tuple(index1.T)].T
            vertex = clus.shape[1] + 1
            dxi = \
                + s01.prod(axis=-1).mean(axis=-1) \
                - s11.prod(axis=-1).mean(axis=-1)
            energy += (dxi * eci) * vertex
        return energy

    def get_exchange_de_small(self, site_s0, site_s1):
        """
        小さな cell 用の計算 (cell のサイズが相関関数のサイズより小さい場合)
        周期境界条件でeciがフリップするサイト自身を参照してしまう場合、
        上の方法 (get_exchange_de) では計算誤差が生じてしまう
        小さな cell の場合 total energy を出すのは早いので、
        単純に、それを使って de を出す
        site1, 2 の spin を交換した場合のエネルギー差を return する
        site1, 2 の表記は x, y, z, site_id 4 要素の array
        """
        de = - self.get_energy()
        self.cell[tuple(site_s0.T)] = 1
        self.cell[tuple(site_s1.T)] = 0
        de += self.get_energy()
        self.cell[tuple(site_s1.T)] = 1
        self.cell[tuple(site_s0.T)] = 0
        return de * self.size ** 3 * 4

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
            np.argwhere(self.cell[..., 0:4] == 0), 's') / self.size

        poss2 = self.conv2real_coordinates(
            np.argwhere(self.cell[..., 0:4] == 1), 's') / self.size

        num1 = poss1.shape[0]

        num2 = poss2.shape[0]

        lines = "mc\n"
        lines += "{0}\n".format(self.size*5)
        lines += "1.00 0 0\n0 1.00 0\n0 0 1.00\n"
        lines += "Cu Al\n"
        lines += "{0} {1}\n".format(num1, num2)
        lines += "Direct\n"
        for site in poss1:
            lines += " ".join([str(x) for x in site]) + "\n"

        for site in poss2:
            lines += " ".join([str(x) for x in site]) + "\n"

        with open(fname, 'w') as wfile:
            wfile.write(lines)

    def make_xdatcar(self, fname):
        poss1 = self.conv2real_coordinates(
            np.argwhere(self.cell[..., 0:4] == 0), 's') / self.size
        # poss2 = self.conv2real_coordinates(
        #     np.argwhere(self.cell[..., 4:8] == 1), 'i') / self.size
        num1 = poss1.shape[0]
        # num2 = poss2.shape[0]
        lines = "\n"
        for site in poss1:
            lines += " ".join([str(x) for x in site]) + "\n"
        # for site in poss2:
        #     lines += " ".join([str(x) for x in site]) + "\n"
        with open(fname, 'a') as wfile:
            wfile.write(lines)

    @staticmethod
    def conv2real_coordinates(coords, s_or_i):
        """
        実空間の座標に変換
        """
        n = {'s': 0, 'i': 4}[s_or_i]
        trans = np.array(NaClSite.SITES)
        real_coords = coords[:, 0:3] + trans[coords[:, 3] + n]
        return real_coords

    def save_cell(self, title):
        """
        cellをpickleデータに保存する
        """
        fname = title + ".pickle"
        with open(fname, 'wb') as wbfile:
            pickle.dump(self.cell, wbfile)

    def from_pickle_cell(self, fname):
        """
        size は揃えておく必要がある
        """
        with open(fname, 'rb') as rbfile:
            self.cell = pickle.load(rbfile)

    def get_conc(self):
        return self.cell.mean()

