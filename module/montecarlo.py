#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
from __future__ import division
import os
import re
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
        """return beta"""
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

    def trans_prob_de2(self, de):
        """
        エネルギー差の行列から遷移確率の行列を返す
        """
        w = np.exp(-de/self.beta)
        return w

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

    def __select_pair_sites2(self, idr=(0,8)):
        """
        s = 0/1 の site の pair をランダムに選択する関数を生成する
        """
        num_s0 = range((self.cell.cell[:, :, :, idr[0]:idr[1]] == 0).sum())
        num_s1 = range((self.cell.cell[:, :, :, idr[0]:idr[1]] == 1).sum())
        def select_pair_sites():
            """
            num_s0, s1 を memo 化しておく
            """
            s0 = self.cell.cell[:, :, :, idr[0]:idr[1]] == 0
            site_s0 = (
                np.array(np.where(s0)))[:, np.random.choice(
                    num_s0, 1, replace=False)[0]]
            s1 = self.cell.cell[:, :, :, idr[0]:idr[1]] == 1
            site_s1 = (
                np.array(np.where(s1)))[:, np.random.choice(
                    num_s1, 1, replace=False)[0]]
            site_s0 += np.array([0, 0, 0, idr[0]])
            site_s1 += np.array([0, 0, 0, idr[0]])
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

    def __loop_bcci_micro_single(self, steps, func):
        return self.__loop_2sub_micro_single(steps, func, (0, 2))

    def __loop_2sub_micro_single(self, steps, func, idr):
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
        select_pair_sites_s = self.__select_pair_sites2(idr=idr)
        select_pair_sites_i = self.__select_pair_sites2(idr=idr)
        # size = self.cell.size**3*4
        delta_e = 0
        size = self.cell.size ** 3 * 8
        # initial = self.cell.get_energy()
        for _ in range(steps):
            # 01s
            site_s0, site_s1 = select_pair_sites_s()
            # 02s
            de = func(site_s0, site_s1)
            p = float(self.trans_prob(de))
            # print(p)
            rand = float(np.random.rand(1))
            s = {False: [0, 1], True: [1, 0]}[rand < p]

            # 03s
            self.cell.cell[tuple(site_s0.T)] = s[0]
            self.cell.cell[tuple(site_s1.T)] = s[1]
            out_de = {False: 0, True: de}[rand < p]
            delta_e += out_de

            # 01
            site_s0, site_s1 = select_pair_sites_i()
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

        return delta_e / size
        # print(delta_e/size)
        # final = self.cell.get_energy()
        # print(final - initial)

    def loop_bcci_micro_single(self, steps):
        """
        """
        return self.__loop_bcci_micro_single(steps, self.cell.get_exchange_de)

    def loop_fcci_micro_single(self, steps):
        return self.__loop_2sub_micro_single(steps, self.cell.get_exchange_de, (0, 4))

    def loop_fcc_micro_single(self, steps):
        """
        normal size の cell 用
        """
        return self.__loop_fcc_micro_single(steps, self.cell.get_exchange_de)

    def loop_fcc_micro_single_small(self, steps):
        """
        fcc cell の spin flip loop
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
        xis = [[self.cell.get_xi_sub(), self.cell.get_xi_int()]]
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
            self.cell.cell[..., 4:8][jmat] = 1 - self.cell.cell[..., 4:8][jmat]

            b = judge.sum()
            energy.append([a, b, self.cell.get_energy()])
            xis.append([self.cell.get_xi_sub(), self.cell.get_xi_int()])
            # print(self.cell.get_energy())
            # print(len(self.cell.get_xi_sub()))
            # print(len(self.cell.get_xi_int()))

        return energy, xis

    def loop_nacl_micro_fix_sub(self, steps):
        """
        NaCl 型 cell の spin flip loop
        micro canonical ensemble
        s-site は固定
        """
        print("initial concentration")
        print(self.cell.get_conc())
        print("initial energy")
        print(self.cell.get_energy())
        print("start")
        energy = [[0, 0, self.cell.get_energy()]]
        xis = [[self.cell.get_xi_sub(), self.cell.get_xi_int()]]
        for _ in range(steps):
            # de = self.cell.get_flip_de_sub().reshape(4, -1).T.reshape(
            #     self.cell.size, self.cell.size, self.cell.size, 4)
            # m2p = self.cell.cell[..., 0:4] == 0
            # p2m = self.cell.cell[..., 0:4] == 1
            # is_m2p_major = m2p.sum() > p2m.sum()
            # major = {True: m2p, False: p2m}[is_m2p_major]
            # minor = {True: p2m, False: m2p}[is_m2p_major]
            # pair = np.random.choice(
            #     range(major.sum()), minor.sum(), replace=False)
            # delta = de[tuple(np.array(np.where(major))[:, pair])] + \
            #     de[tuple(np.array(np.where(minor)))]
            # p = self.trans_prob(delta)
            # rand = np.random.rand(p.size)
            # judge = (rand < p)
            # jmat = major * False
            # jmat[tuple(np.array(np.where(major))[:, pair])] = judge
            # jmat[tuple(np.array(np.where(minor)))] = judge
            # self.cell.cell[..., 0:4][jmat] = 1 - self.cell.cell[..., 0:4][jmat]

            a = 0

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
            self.cell.cell[..., 4:8][jmat] = 1 - self.cell.cell[..., 4:8][jmat]

            b = judge.sum()
            energy.append([a, b, self.cell.get_energy()])
            xis.append([self.cell.get_xi_sub(), self.cell.get_xi_int()])
            # print(self.cell.get_energy())
            # print(len(self.cell.get_xi_sub()))
            # print(len(self.cell.get_xi_int()))

        return energy, xis

    def loop_nacl_micro_fix_int(self, steps):
        """
        NaCl 型 cell の spin flip loop
        micro canonical ensemble
        s-site は固定
        """
        print("initial concentration")
        print(self.cell.get_conc())
        print("initial energy")
        print(self.cell.get_energy())
        print("start")
        energy = [[0, 0, self.cell.get_energy()]]
        xis = [[self.cell.get_xi_sub(), self.cell.get_xi_int()]]
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

            # de = self.cell.get_flip_de_int().reshape(4, -1).T.reshape(
            #     self.cell.size, self.cell.size, self.cell.size, 4)
            # m2p = self.cell.cell[..., 4:8] == 0
            # p2m = self.cell.cell[..., 4:8] == 1
            # is_m2p_major = m2p.sum() > p2m.sum()
            # major = {True: m2p, False: p2m}[is_m2p_major]
            # minor = {True: p2m, False: m2p}[is_m2p_major]
            # pair = np.random.choice(
            #     range(major.sum()), minor.sum(), replace=False)
            # delta = de[tuple(np.array(np.where(major))[:, pair])] + \
            #     de[tuple(np.array(np.where(minor)))]
            # p = self.trans_prob(delta)
            # rand = np.random.rand(p.size)
            # judge = (rand < p)
            # jmat = major * False
            # jmat[tuple(np.array(np.where(major))[:, pair])] = judge
            # jmat[tuple(np.array(np.where(minor)))] = judge
            # self.cell.cell[..., 4:8][jmat] = 1 - self.cell.cell[..., 4:8][jmat]

            b = 0
            energy.append([a, b, self.cell.get_energy()])
            xis.append([self.cell.get_xi_sub(), self.cell.get_xi_int()])
            # print(self.cell.get_energy())
            # print(len(self.cell.get_xi_sub()))
            # print(len(self.cell.get_xi_int()))

        return energy, xis

    def loop_bcci_micro(self, steps):
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
        # xis = [[self.cell.get_xi_sub(), self.cell.get_xi_int()]]
        for _ in range(steps):
            de = self.cell.get_flip_de()[:,:,:,:2].reshape(2, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 2)
            m2p = self.cell.cell[..., 0:2] == 0
            p2m = self.cell.cell[..., 0:2] == 1
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
            self.cell.cell[..., 0:2][jmat] = 1 - self.cell.cell[..., 0:2][jmat]

            a = judge.sum()

            de = self.cell.get_flip_de()[:,:,:,2:].reshape(6, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 6)
            m2p = self.cell.cell[..., 2:8] == 0
            p2m = self.cell.cell[..., 2:8] == 1
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
            self.cell.cell[..., 2:8][jmat] = 1 - self.cell.cell[..., 2:8][jmat]

            b = judge.sum()
            energy.append([a, b, self.cell.get_energy()])
            # print(self.cell.get_energy())
            # print(len(self.cell.get_xi_sub()))
            # print(len(self.cell.get_xi_int()))
        return energy

    def loop_bcci_micro_fix_sub(self, steps):
        """
        NaCl 型 cell の spin flip loop
        micro canonical ensemble
        s-site は固定
        """
        print("initial concentration")
        print(self.cell.get_conc())
        print("initial energy")
        print(self.cell.get_energy())
        print("start")
        energy = [[0, 0, self.cell.get_energy()]]
        xis = [[self.cell.get_xi_sub(), self.cell.get_xi_int()]]
        for _ in range(steps):
            # de = self.cell.get_flip_de_sub().reshape(4, -1).T.reshape(
            #     self.cell.size, self.cell.size, self.cell.size, 4)
            # m2p = self.cell.cell[..., 0:4] == 0
            # p2m = self.cell.cell[..., 0:4] == 1
            # is_m2p_major = m2p.sum() > p2m.sum()
            # major = {True: m2p, False: p2m}[is_m2p_major]
            # minor = {True: p2m, False: m2p}[is_m2p_major]
            # pair = np.random.choice(
            #     range(major.sum()), minor.sum(), replace=False)
            # delta = de[tuple(np.array(np.where(major))[:, pair])] + \
            #     de[tuple(np.array(np.where(minor)))]
            # p = self.trans_prob(delta)
            # rand = np.random.rand(p.size)
            # judge = (rand < p)
            # jmat = major * False
            # jmat[tuple(np.array(np.where(major))[:, pair])] = judge
            # jmat[tuple(np.array(np.where(minor)))] = judge
            # self.cell.cell[..., 0:4][jmat] = 1 - self.cell.cell[..., 0:4][jmat]

            a = 0

            de = self.cell.get_flip_de_int().reshape(4, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 4)
            m2p = self.cell.cell[..., 2:8] == 0
            p2m = self.cell.cell[..., 2:8] == 1
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
            self.cell.cell[..., 2:8][jmat] = 1 - self.cell.cell[..., 2:8][jmat]

            b = judge.sum()
            energy.append([a, b, self.cell.get_energy()])
            xis.append([self.cell.get_xi_sub(), self.cell.get_xi_int()])
            # print(self.cell.get_energy())
            # print(len(self.cell.get_xi_sub()))
            # print(len(self.cell.get_xi_int()))
        return energy, xis

    def loop_bcci_micro_fix_int(self, steps):
        """
        NaCl 型 cell の spin flip loop
        micro canonical ensemble
        s-site は固定
        """
        print("initial concentration")
        print(self.cell.get_conc())
        print("initial energy")
        print(self.cell.get_energy())
        print("start")
        energy = [[0, 0, self.cell.get_energy()]]
        xis = [[self.cell.get_xi_sub(), self.cell.get_xi_int()]]
        for _ in range(steps):
            de = self.cell.get_flip_de_sub().reshape(4, -1).T.reshape(
                self.cell.size, self.cell.size, self.cell.size, 4)
            m2p = self.cell.cell[..., 0:2] == 0
            p2m = self.cell.cell[..., 0:2] == 1
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
            self.cell.cell[..., 0:2][jmat] = 1 - self.cell.cell[..., 0:2][jmat]

            a = judge.sum()

            # de = self.cell.get_flip_de_int().reshape(4, -1).T.reshape(
            #     self.cell.size, self.cell.size, self.cell.size, 4)
            # m2p = self.cell.cell[..., 4:8] == 0
            # p2m = self.cell.cell[..., 4:8] == 1
            # is_m2p_major = m2p.sum() > p2m.sum()
            # major = {True: m2p, False: p2m}[is_m2p_major]
            # minor = {True: p2m, False: m2p}[is_m2p_major]
            # pair = np.random.choice(
            #     range(major.sum()), minor.sum(), replace=False)
            # delta = de[tuple(np.array(np.where(major))[:, pair])] + \
            #     de[tuple(np.array(np.where(minor)))]
            # p = self.trans_prob(delta)
            # rand = np.random.rand(p.size)
            # judge = (rand < p)
            # jmat = major * False
            # jmat[tuple(np.array(np.where(major))[:, pair])] = judge
            # jmat[tuple(np.array(np.where(minor)))] = judge
            # self.cell.cell[..., 4:8][jmat] = 1 - self.cell.cell[..., 4:8][jmat]

            b = 0
            energy.append([a, b, self.cell.get_energy()])
            xis.append([self.cell.get_xi_sub(), self.cell.get_xi_int()])
            # print(self.cell.get_energy())
            # print(len(self.cell.get_xi_sub()))
            # print(len(self.cell.get_xi_int()))

        return energy, xis


class BcciOctaSite(object):
    """
    Bcc + Octahedron-interstital 8 点のサイト情報を記録する
    s-site 2 点 [0:2]
    i-site 6 点 [3:8]
    args:
        spins: spin の配列を記した 8 要素の list
    cls variables:
        SPIN: spin 軸の取り方 default: [1, 0]
        SITES: 各サイトの座標
        CONC_SUB: s-site の SPIN[0] の濃度
        CONC_INT: i-site の SPIN[0] の濃度
        EQUIV: 対称操作が等価な site
    """
    SPIN = [1, 0]
    _site_s0 = [0.0, 0.0, 0.0]
    _site_s1 = [1/2, 1/2, 1/2]
    _site_i0 = [1/2, 0.0, 0.0]
    _site_i1 = [0.0, 1/2, 0.0]
    _site_i2 = [0.0, 0.0, 1/2]
    _site_i3 = [0.0, 1/2, 1/2]
    _site_i4 = [1/2, 0.0, 1/2]
    _site_i5 = [1/2, 1/2, 0.0]
    SITES = [_site_s0, _site_s1, _site_i0, _site_i1,
             _site_i2, _site_i3, _site_i4, _site_i5]
    CONC_INT = 0.5
    CONC_SUB = 0.5
    LABEL = {0: 'C', 1: 'C', 2: 'A', 3: 'A',
             4: 'A', 5: 'A', 6: 'A', 7: 'A'}
    EQUIV = {0: 'C', 1: 'C', 2: 'A1', 3: 'A2',
             4: 'A3', 5: 'A1', 6: 'A2', 7: 'A3'}
    COUNTER = [2, 2, 6, 6, 6, 6, 6, 6]

    def __init__(self, spins):
        self.spins = spins

    @classmethod
    def conv_site2idx(cls, site):
        """
        座標の直接表記を class 内の index 表記に変換して return する
        """
        index = [math.floor(x) for x in site]
        decim = [x - math.floor(x) for x in site]
        index.append({tuple(cls.SITES[i]): i
                      for i in range(len(cls.SITES))}[tuple(decim)])
        return index

    @classmethod
    def random(cls, conc):
        """
        spin を site に random に配置する
        cls.CONC_INT と cls.CONC_SUB で i-, s-site の濃度を指定する
        """
        if conc:
            cls.CONC_INT = conc[0]
            cls.CONC_SUB = conc[1]
        i = cls.CONC_INT
        s = cls.CONC_SUB
        spins = [cls.SPIN[int(round(random.uniform(0.5-s/2., 1-s/2.)))]
                 for _ in range(2)] + \
                [cls.SPIN[int(round(random.uniform(0.5-i/2., 1-i/2.)))]
                 for _ in range(6)]
        return cls(spins)

    @classmethod
    def BCC_u(cls, _):
        """
        BCC 構造
        s-site が up-spin
        """
        i, j = cls.SPIN
        spins = [i, i, j, j, j, j, j, j]
        return cls(spins)

    @classmethod
    def BCC_d(cls, _):
        """
        BCC 構造
        s-site が up-spin
        """
        i, j = cls.SPIN
        spins = [j, j, j, j, j, j, j, j]
        return cls(spins)

    @classmethod
    def A3C(cls, _):
        i, j = cls.SPIN
        spins = [i, i, i, i, i, i, i, i]
        return cls(spins)

    @classmethod
    def A3D(cls, _):
        i, j = cls.SPIN
        spins = [j, j, i, i, i, i, i, i]
        return cls(spins)

    @classmethod
    def A2D(cls, _):
        i, j = cls.SPIN
        spins = [j, j, i, j, i, i, j, i]
        return cls(spins)

    @classmethod
    def A1D(cls, _):
        i, j = cls.SPIN
        spins = [j, j, j, i, j, j, i, j]
        return cls(spins)

    @classmethod
    def A2C(cls, _):
        i, j = cls.SPIN
        spins = [i, i, j, i, i, j, i, i]
        return cls(spins)

    @classmethod
    def A1C(cls, _):
        i, j = cls.SPIN
        spins = [i, i, j, j, i, j, j, i]
        return cls(spins)

    @classmethod
    def ACD(cls, _):
        i, j = cls.SPIN
        spins = [j, i, i, j, j, i, j, j]
        return cls(spins)


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
    def conv_site2idx(cls, site):
        """
        座標の直接表記を class 内の index 表記に変換して return する
        """
        index = [math.floor(x) for x in site]
        decim = [x - math.floor(x) for x in site]
        index.append({tuple(cls.SITES[i]): i
                      for i in range(len(cls.SITES))}[tuple(decim)])
        return index

    @classmethod
    def random(cls, _):
        """
        spin を site に random に配置する
        cls.CONC_INT と cls.CONC_SUB で i-, s-site の濃度を指定する
        """
        i = cls.CONC_INT
        s = cls.CONC_SUB
        spins = [cls.SPIN[int(round(random.uniform(0.5-s/2., 1-s/2.)))]
                 for _ in range(4)] + \
                [cls.SPIN[int(round(random.uniform(0.5-i/2., 1-i/2.)))]
                 for _ in range(4)]
        return cls(spins)

    @classmethod
    def FCC_u(cls, _):
        """
        FCC 構造
        s-site が up-spin
        """
        i, j = cls.SPIN
        spins = [i, i, i, i, j, j, j, j]
        return cls(spins)

    @classmethod
    def FCC_d(cls, _):
        """
        FCC 構造
        s-site が down-spin
        """
        _, j = cls.SPIN
        spins = [j, j, j, j, j, j, j, j]
        return cls(spins)

    @classmethod
    def L10(cls, _):
        """
        L10 構造
        """
        i, j = cls.SPIN
        spins = [i, i, j, j, j, j, j, j]
        return cls(spins)

    @classmethod
    def L12_A(cls, _):
        """
        A3B L12 構造
        """
        i, j = cls.SPIN
        spins = [i, i, i, j, j, j, j, j]
        return cls(spins)

    @classmethod
    def NaCl_uu(cls, _):
        """
        NaCl 型構造
        i-, s-site ともに up-spin
        """
        i, _ = cls.SPIN
        spins = [i, i, i, i, i, i, i, i]
        return cls(spins)

    @classmethod
    def NaCl_du(cls, _):
        """
        NaCl型構造
        i-site down-spin, s-site up-spin
        """
        i, j = cls.SPIN
        spins = [j, j, j, j, i, i, i, i]
        return cls(spins)

    @classmethod
    def order_00(cls, _):
        """
        NaCl型構造
        i-site down-spin, s-site up-spin
        """
        i, j = cls.SPIN
        spins = [j, j, j, i, j, j, i, j]
        return cls(spins)

    @classmethod
    def order_01(cls, _):
        """
        NaCl型構造
        i-site down-spin, s-site up-spin
        """
        i, j = cls.SPIN
        spins = [i, i, i, j, j, j, i, j]
        return cls(spins)

    def __str__(self):
        return str(self.spins)

    def __getitem__(self, index):
        return self.spins[index]


class CellIOMixin(object):
    """
    cell の読み書き用の mixin
    """
    def save_cell(self, title):
        """
        cellをpickleデータに保存する
        """
        fname = title + ".pickle"
        with open(fname, 'wb') as wbfile:
            pickle.dump(self.cell, wbfile)

    def load_cell(self, fname):
        """
        size は揃えておく必要がある
        """
        with open(fname, 'rb') as rbfile:
            self.cell = pickle.load(rbfile)


class NaClXtal(CellIOMixin):
    """
    NaCl 型の結晶格子を作成する
    各サイトは index 表記になっているため、座標の参照の仕方が面倒
    4 番目のコラムが unit cell 内の site の index

    0 以外の site から参照する際に self.TRANS を利用する
    例として、site-1 から site-2 を参照させる場合
    self.TRANS[1, 2] を index に足し合わせる
    (直感的に分かりづらいのが問題)

    args:
        size: cell の size (total の 原子数は size ** 3 × 8sites)
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
    def __init__(self, size, arrange='random', conc=(0.5, 0.5), ecis=None):
        self.size = size
        NaClSite.CONC_SUB = conc[0]
        NaClSite.CONC_INT = conc[1]
        mode = {'random': NaClSite.random, 'L10': NaClSite.L10,
                'FCC_u': NaClSite.FCC_u, 'NaCl_uu': NaClSite.NaCl_uu,
                'FCC_d': NaClSite.FCC_d, 'NaCl_du': NaClSite.NaCl_du,
                'order_00': NaClSite.order_00, 'order_01': NaClSite.order_01,
               }[arrange]
        self.cell = np.array([[[mode(conc).spins for z in range(size)]
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

    @staticmethod
    def get_trans():
        """
        相対座標へ変換するための matrix を return する
        """
        sites = np.array(NaClSite.SITES)
        size = sites.shape[0]
        tmp = sites.reshape(size, 1, 3) + sites.reshape(1, size, 3)
        index = np.array([[NaClSite.conv_site2idx(x) for x in y] for y in tmp])
        obj = np.c_[
            np.zeros_like(np.arange(size*3).reshape(size, 3)), np.arange(size)]
        return index - obj

    @classmethod
    def from_pickle_ecis(cls, fname, size=10,
                         arrange='random', conc=(0.5, 0.5)):
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
        return cls(size, arrange, conc, ecis_data)

    def get_xi_sub(self):
        """
        置換型サイトの相関関数を出力
        """
        xis = []
        for idxs in self.sub_idxs:
            spins = self.cell[idxs].T
            xi = spins.prod(axis=-1).mean(axis=-1)
            xis.append(xi.mean())
        return xis

    def get_xi_int(self):
        """
        侵入型サイトの相関関数を出力
        """
        xis = []
        for idxs in self.int_idxs:
            spins = self.cell[idxs].T
            xi = spins.prod(axis=-1).mean(axis=-1)
            xis.append(xi.mean())
        return xis

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
            # print(len(idxs))
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

    # def get_exchange_de(self, site_s0, site_s1):
    #     """
    #     spin 交換によるエネルギー差
    #     site_s0 は s = 0 のサイト
    #     site_s1 は s = 1 のサイトを指定
    #     site_s0, s1 の spin を交換した場合のエネルギー差を return する
    #     site_s0, s1 の表記は [x, y, z, site_id] の 4 要素の array
    #     othrer variables:
    #         tmp_cell: spin_flip 後の cell
    #     """
    #     tmp_cell = self.cell.copy()
    #     tmp_cell[tuple(site_s0.T)] = 1
    #     tmp_cell[tuple(site_s1.T)] = 0

    #     exchange_de = 0
    #     for clus, eci in zip(self.clusters, self.ecis):
    #         index1 = (clus + self.TRANS[site_s0[3]][clus[..., 3]] +
    #                   np.r_[site_s0[0:3], [0]])
    #         index1[..., 0:3] %= self.size
    #         index2 = (clus + self.TRANS[site_s1[3]][clus[..., 3]] +
    #                   np.r_[site_s1[0:3], [0]])
    #         index2[..., 0:3] %= self.size

    #         # s=1 のサイトの相関関数
    #         s11 = self.cell[tuple(index2.T)].T
    #         # s=0 だったところに 1 が入った時の相関関数
    #         # この時、 s=1 だったところは 0になっているので
    #         # temp_cell を使う
    #         s01 = tmp_cell[tuple(index1.T)].T
    #         vertex = clus.shape[1] + 1
    #         dxi = \
    #             + s01.prod(axis=-1).mean(axis=-1) \
    #             - s11.prod(axis=-1).mean(axis=-1)
    #         exchange_de += (dxi * eci) * vertex
    #     return exchange_de

    def get_conc(self):
        """
        組成比を return
        return:
            {'inter': , 'subs': }
        """
        return {'subs': self.cell[..., 0:4].mean(),
                'inter': self.cell[..., 4:8].mean()}

    def make_poscar(self, fname, s_spin, i_spin, with_s_oposit=False):
        """
        self.cellをPOSCARのformatで出力する
        attributes
            fname: 出力のファイル名
            s_spin: 書き出す substitutional site の spin の値
            i_spin: 書き出す interstitial site の spin の値
        other variables
            poss1, poss2: 元素A, Bの座標
            num1, num2: 元素A, Bのtotalの元素数
            lines: 出力
        """
        poss1 = self.conv2real_coordinates(
            np.argwhere(self.cell[..., 0:4] == s_spin), 's') / self.size
        poss2 = self.conv2real_coordinates(
            np.argwhere(self.cell[..., 4:8] == i_spin), 'i') / self.size
        num1 = poss1.shape[0]
        num2 = poss2.shape[0]
        poss_o = None
        num_o = 0
        if with_s_oposit:
            poss_o = self.conv2real_coordinates(
                np.argwhere(self.cell[..., 0:4] == 1 - s_spin), 's') / self.size
            num_o = poss_o.shape[0]
        lines = "mc\n"
        lines += "{0}\n".format(self.size*2.5)
        lines += "1.00 0 0\n0 1.00 0\n0 0 1.00\n"
        lines += "Cr C Fe\n"
        lines += "{0} {1} {2}\n".format(num1, num2, num_o)
        lines += "Direct\n"
        for site in poss1:
            lines += " ".join([str(x) for x in site]) + "\n"
        for site in poss2:
            lines += " ".join([str(x) for x in site]) + "\n"
        if with_s_oposit:
            for site in poss_o:
                lines += " ".join([str(x) for x in site]) + "\n"
        with open(fname, 'w') as wfile:
            wfile.write(lines)

    def load_poscar(self, fname, s_spin, i_spin):
        """
        書き出した poscar を読み込む
        size は object の self.size と合致していなければならない
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        if self.size != int(float(lines[1]) / 2.5):
            print("size is different !")
            return
        num_s, num_i, _ = (int(x) for x in lines[6].split())
        for s in range(num_s):
            site = [float(x) * self.size for x in lines[s+8].split()]
            site_idx = np.array(self.conv_NaCl(site))
            self.cell[tuple(site_idx.T)] = s_spin
        for i in range(num_i):
            site = [float(x) * self.size for x in lines[i+8+num_s].split()]
            site_idx = np.array(self.conv_NaCl(site))
            self.cell[tuple(site_idx.T)] = i_spin

    @staticmethod
    def conv2real_coordinates(coords, s_or_i):
        """
        実空間の座標に変換
        """
        n = {'s': 0, 'i': 4}[s_or_i]
        trans = np.array(NaClSite.SITES)
        real_coords = coords[:, 0:3] + trans[coords[:, 3] + n]
        return real_coords

    @staticmethod
    def conv_NaCl(site):
        """
        NaCl site の表記に変換する
        dcimal to integer
        """
        integ = []
        decim = []
        for coord in site:
            tmp = math.floor(round(coord*2, 1) / 2)
            integ.append(tmp)
            decim.append(round((coord - tmp), 2))
        site_id = {(0.0, 0.0, 0.0): 0, (0.0, 0.5, 0.5): 1,
                   (0.5, 0.0, 0.5): 2, (0.5, 0.5, 0.0): 3,
                   (0.5, 0.0, 0.0): 4, (0.0, 0.5, 0.0): 5,
                   (0.0, 0.0, 0.5): 6, (0.5, 0.5, 0.5): 7}[tuple(decim)]
        integ.append(site_id)
        return integ

    def get_exchange_de(self, site_s0, site_s1):
        size = self.size ** 4 * 8
        e0 = self.get_energy()
        self.cell[tuple(site_s0.T)] = 1
        self.cell[tuple(site_s1.T)] = 0
        e1 = self.get_energy()
        self.cell[tuple(site_s0.T)] = 0
        self.cell[tuple(site_s1.T)] = 1
        return (e1 - e0) * size


class SubLattXtal(object):
    """
    二副格子型の結晶格子を作成する
    各サイトは index 表記になっているため、座標の参照の仕方が面倒
    4 番目のコラムが unit cell 内の site の index

    0 以外の site から参照する際に self.TRANS を利用する
    例として、site = 1 から site = 2 を参照させる場合
    self.TRANS[1, 2] を index に足し合わせる
    (直感的に分かりづらいのが問題)

    args:
        size: cell のsize size^3*8sites を作成する
        arrange: 原子配列の初期値
                 選択肢は mode を参照
        conc: up-spinの組成比の初期値
              arrange=random のみ使用
              (第一副格子, 第二副格子) の順
        ecis: ECIs
              format は list形式
              第一、第二副格子の順で格納した
              dict形式eciの情報を持つ (CEMParserから作成する)
    variables:
        nodup_int_ecis: i-site の ecis で s-site と重複していないもののlist
                        total energy を算出する際に使用
                        これをやらないと重複して相関関数を数えてしまうので、
                        total energy が変わってしまう
        TRANS: 相対座標
               原点を別の座標にとった際、他の座標を参照する際に使用
               [新しい原点の座標, 参照する座標]として記述している
               例: site3からsite2を参照する場合
                   + [1, 0, 0, -1]した座標がそのサイトになる
        COORDS: 全ての格子点の座標 (site_id は 0 のみ)
    """
    def __init__(self, size, ecis, site_class, arrange='random', conc=None):
        self.size = size
        mode = {'random': site_class.random,
                'BCC_u': site_class.BCC_u,
                'BCC_d': site_class.BCC_d,
                'A3C': site_class.A3C, 'A3D': site_class.A3D,
                'A2C': site_class.A2C, 'A2D': site_class.A2D,
                'A1C': site_class.A1C, 'A1D': site_class.A1D,
                'ACD': site_class.ACD}[arrange]

        self.cell = np.array([[[mode(conc).spins for z in range(size)]
                               for y in range(size)] for x in range(size)])
        self.site_class = site_class
        self.COORDS = np.fromiter(chain.from_iterable(
            product(range(0, self.size), repeat=3)), np.int).reshape(-1, 3)
        self.COORDS = np.c_[self.COORDS,
                            np.zeros_like(
                                self.COORDS[:, 0], np.int)].reshape(-1, 1, 1, 4)
        self.TRANS = self.get_trans()
        # print(self.TRANS)

        self.clusters = [ecis[site_class.EQUIV[i]][0][1:]
                         for i in range(len(site_class.SITES))]
        self.ecis = [ecis[site_class.EQUIV[i]][1][1:]
                     for i in range(len(site_class.SITES))]
        self.eci_point = [ecis[site_class.EQUIV[i]][1][0]
                          for i in range(len(site_class.SITES))]
        # self.idxs = [self.get_clus_idxs(clusters)
        #              for clusters in self.clusters]
        self.vertex = []
        self.num_sublatt = []
        self.sublatt1s = [x for x in site_class.LABEL
                     if site_class.LABEL[x] == site_class.LABEL[0]]
        for i in site_class.LABEL:
            label = [self.TRANS[i][x][3] + x for x in self.sublatt1s]
            tmp = []
            tmp2 = []
            for clus in self.clusters[i]:
                num = 1
                if i in self.sublatt1s:
                    num += len([x for x in clus[0] if x[3] in label])
                else:
                    num += len(clus[0]) - len([x for x in clus[0]
                                               if x[3] in label])
                tmp.append(num)
                tmp2.append(len(clus[0]) + 1)
            self.num_sublatt.append(tmp)
            self.vertex.append(tmp2)

    @classmethod
    def from_pickle_ecis(
        cls, size, fname, site_class, arrange='random', conc=None):
        """
        ecis や cluster の情報が入った pickle data から cls を生成
        load する pickle data の format は
        dict{clusters: [(n - 2) × degen × (vertex - 1) × 4 座標]
             ecis: [(n - 2) * ecis]
             eci_point: point cluster のみの ecis}
        n - 2 は null と point を除いてあることに
        vertex-1 は 原点 [0, 0, 0, 0] を除いてあることに対応する
        """
        with open(fname, 'rb') as rbfile:
            ecis_data = pickle.load(rbfile)
        return cls(size, ecis_data, site_class, arrange, conc)


    def get_exchange_de(self, site_s0, site_s1):
        """
        cell 全体のエネルギーを使って
        spin 交換に伴うエネルギー差を計算
        計算コストでかい・・・
        """
        size = self.size ** 4 * 8
        e0 = self.get_energy()
        self.cell[tuple(site_s0.T)] = 1
        self.cell[tuple(site_s1.T)] = 0
        e1 = self.get_energy()
        self.cell[tuple(site_s0.T)] = 0
        self.cell[tuple(site_s1.T)] = 1
        return (e1 - e0) * size

    def get_trans(self):
        """
        相対座標へ変換するための matrix を return する
        """
        sites = np.array(self.site_class.SITES)
        size = sites.shape[0]
        tmp = sites.reshape(size, 1, 3) + sites.reshape(1, size, 3)
        index = np.array([[self.site_class.conv_site2idx(x)
                           for x in y] for y in tmp])
        obj = np.c_[
            np.zeros_like(np.arange(size*3).reshape(size, 3)), np.arange(size)]
        return index - obj

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
            for site_id in range(8):
                pos = clus + self.TRANS[site_id][clus[..., 3]] + \
                    self.COORDS
                tmp_pos.append(pos)
            index = np.array(tmp_pos).reshape(
                8*self.size**3, size[0], size[1], 4)
            index[..., 0:3] %= self.size
            indexes.append(tuple(index.T))
        return indexes

    # def get_energy(self):
    #     """
    #     siteに仮想的に spin=1 を入れて、その周囲の相関関数を計算
    #     ecis との積から サイト毎のエネルギーを計算する
    #     """
    #     energy = np.zeros_like(np.arange(self.size**3*8))
    #     dummy = np.zeros_like(
    #         np.arange(self.size**3*8).
    #         reshape(8, self.size, self.size, self.size))
    #     for i in range(len(self.ecis)):
    #         for idxs, eci, num_sub, vertex in zip(self.idxs[i],
    #                                               self.ecis[i],
    #                                               self.num_sublatt[i],
    #                                               self.vertex[i]):
    #             correction = num_sub / vertex
    #             spins = self.cell[idxs].T
    #             xi = spins.prod(axis=-1).mean(axis=-1)
    #             dummy00 = dummy.copy()
    #             dummy00[i, :, :, :] = 1
    #             dummy00 = dummy00.reshape(self.size**3*8)
    #             energy = (xi * eci * dummy00) * correction + energy
    #     energy = energy.reshape(8, -1).T.reshape(
    #         self.size, self.size, self.size, 8)
    #     energy_sum = 0
    #     for i in range(len(self.ecis)):
    #         energy[:, :, :, i] += self.eci_point[i]
    #         energy_sum += (energy * self.cell)[:, :, :, i].mean() / \
    #             self.site_class.COUNTER[i]
    #     return energy_sum

    def get_exchange_de2(self, site_s0, site_s1):
        """
        cell 全体のエネルギーを使って
        spin 交換に伴うエネルギー差を計算
        """
        i = site_s0[3]
        j = site_s1[3]
        tmp_cell = self.cell.copy()
        tmp_cell[tuple(site_s0.T)] = 1
        tmp_cell[tuple(site_s1.T)] = 0
        e0 = 0
        e1 = 0
        correction0 = len(self.site_class.SITES) / self.site_class.COUNTER[i]
        for clus, eci, num_sub, vertex in zip(self.clusters[i],
                                              self.ecis[i],
                                              self.num_sublatt[i],
                                              self.vertex[i]):
            correction = num_sub / vertex * correction0
            index0 = (clus + self.TRANS[site_s0[3]][clus[..., 3]] +
                      np.r_[site_s0[0:3], [0]])
            index0[..., 0:3] %= self.size
            # s=0 だったところに 1 が入った時の相関関数
            xi0 = tmp_cell[tuple(index0.T)].T.prod(axis=-1).mean(axis=-1)
            e0 += (xi0 * eci) * vertex * correction
        for clus, eci, num_sub, vertex in zip(self.clusters[j],
                                              self.ecis[j],
                                              self.num_sublatt[j],
                                              self.vertex[j]):
            correction = num_sub / vertex * correction0
            index1 = (clus + self.TRANS[site_s1[3]][clus[..., 3]] +
                      np.r_[site_s1[0:3], [0]])
            index1[..., 0:3] %= self.size
            # s=1 のサイトの相関関数
            xi1 = self.cell[tuple(index1.T)].T.prod(axis=-1).mean(axis=-1)
            e1 += (xi1 * eci) * vertex * correction
        return e0 - e1

    def get_energy(self):
        """
        s = 1 のところのみを抽出して計算する
        希薄濃度域ではこちらが速いと思われる

        希薄濃度域でなくてもこの方法が速い
        過去の計算方式を刷新する
        """
        s1 = self.cell[:, :, :, :] == 1
        site_s1 = (np.array(np.where(s1)))
        site_s1_i0 = (np.array(np.where(site_s1[3]==0)))
        e1 = 0.0
        for i in set(site_s1[3]):
            correction0 = len(self.site_class.SITES) / self.site_class.COUNTER[i]
            site_s1_i = (np.array(np.where(site_s1[3] == i)))
            e1 += self.eci_point[i] * len(site_s1_i[0]) * correction0
            for clus, eci, num_sub, vertex in zip(self.clusters[i],
                                                  self.ecis[i],
                                                  self.num_sublatt[i],
                                                  self.vertex[i]):
                correction = num_sub / vertex * correction0
                index1 = ((clus + self.TRANS[i][clus[..., 3]]).reshape(1, -1, vertex-1, 4) +
                          np.r_[site_s1[:, tuple(site_s1_i)][0:3].reshape(3, -1), np.array([0]*len(site_s1_i[0])).reshape(1, -1)].T.reshape(-1,1,1,4))
                index1[..., 0:3] %= self.size
                xi1 = self.cell[tuple(index1.T)].T.prod(axis=-1).mean(axis=-1).sum(axis=-1)
                e1 += (xi1 * eci) * correction
        return e1 / self.size**3 / 8

    def get_flip_de(self):
        """
        置換型サイトの spin flip によるエネルギー差をサイト毎に return
        """
        de = 0
        for idxs, eci in zip(self.idxs, self.ecis):
            print(len(idxs))
            print("pass")
            spins = self.cell[idxs].T
            xi = spins.prod(axis=-1).mean(axis=-1)
            spins_inv = spins
            spins_inv[..., 0] = 1 - spins[..., 0]
            xi_inv = spins_inv.prod(axis=-1).mean(axis=-1)
            dxi = xi_inv - xi
            de = dxi * eci + de
        return de

    def get_conc(self):
        """
        組成比を return
        return:
            {'inter': , 'subs': }
        """
        return {'subs': self.cell[..., 0:2].mean(),
                'inter': self.cell[..., 2:8].mean()}

    @staticmethod
    def conv2real_coordinates(coords, s_or_i):
        """
        実空間の座標に変換
        """
        n = {'s': 0, 'i': 2}[s_or_i]
        trans = np.array(BcciOctaSite.SITES)
        real_coords = coords[:, 0:3] + trans[coords[:, 3] + n]
        return real_coords

    def make_poscar(self, fname, s_spin, i_spin, with_s_oposit=False):
        """
        self.cellをPOSCARのformatで出力する
        attributes
            fname: 出力のファイル名
            s_spin: 書き出す substitutional site の spin の値
            i_spin: 書き出す interstitial site の spin の値
        other variables
            poss1, poss2: 元素A, Bの座標
            num1, num2: 元素A, Bのtotalの元素数
            lines: 出力
        """
        poss1 = self.conv2real_coordinates(
            np.argwhere(self.cell[..., 0:2] == s_spin), 's') / self.size
        poss2 = self.conv2real_coordinates(
            np.argwhere(self.cell[..., 2:8] == i_spin), 'i') / self.size
        num1 = poss1.shape[0]
        num2 = poss2.shape[0]
        poss_o = None
        num_o = 0
        if with_s_oposit:
            poss_o = self.conv2real_coordinates(
                np.argwhere(self.cell[..., 0:4] == 1 - s_spin), 's') / self.size
            num_o = poss_o.shape[0]
        lines = "mc\n"
        lines += "{0}\n".format(self.size*2.5)
        lines += "1.00 0 0\n0 1.00 0\n0 0 1.00\n"
        lines += "Cr C Fe\n"
        lines += "{0} {1} {2}\n".format(num1, num2, num_o)
        lines += "Direct\n"
        for site in poss1:
            lines += " ".join([str(x) for x in site]) + "\n"
        for site in poss2:
            lines += " ".join([str(x) for x in site]) + "\n"
        if with_s_oposit:
            for site in poss_o:
                lines += " ".join([str(x) for x in site]) + "\n"
        with open(fname, 'w') as wfile:
            wfile.write(lines)


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
    def conv_site2idx(cls, site):
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
    np.array object  3 dim × 4 sites
    size: cellの大きさ (total の原子数 size**3 × 4 sites)
    arrange: 'random', 'L12', 'L10' などを指定
    conc: A, B 二元系における A の濃度の初期値 (random以外では不使用)
    other variables
        TRANS: 相対座標
               0 以外の site から参照する際に self.TRANS を利用する
               例として、site = 1 から site = 2 を参照させる場合
               self.TRANS[1, 2] を index に足し合わせる
               (直感的に分かりづらいのが問題)

        COORDS: 全ての格子点の座標 (site_id は 0 のみ)

    相関関数の計算方法は利便性のために原点 [0, 0, 0, 0] には
    必ず spin = 1 を置いて計算する (xi_o とする)
    各サイトの xi_o を予め計算しておくことで、
    実際のサイト上の xi_o は サイト上の spin * xi_o で計算できる
    spin flip のエネルギー算出の際に計算を簡略化できる

    xi_o 算出のために point cluster は initialize の時点で分離
    (get_clus.py の時点で分離しておく)
    また 各クラスターの原点の座標 [0, 0, 0, 0] も取り除く

    上のやり方は spin 交換での判定の場合は誤差が出てくるので使えない
    (grand canonical では利用できるか)
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
        n - 2 は null と point を除いてあることに
        vertex-1 は 原点 [0, 0, 0, 0] を除いてあることに対応する
        """
        with open(fname, 'rb') as rbfile:
            ecis_data = pickle.load(rbfile)
        return cls(size, arrange, conc, ecis_data)

    def get_clus_idxs(self, clusters):
        """
        fancy index に用いる cluster の index を list 化して return する
        index の shape は
        cluster 種 × 参照する座標 × クラスターの頂点 ×
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
        index = np.array([[QuadSite.conv_site2idx(x) for x in y] for y in tmp])
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
        siteに仮想的に spin = 1 を入れて、その周囲の相関関数を計算
        また、頂点の数との積と取っている (全体のエネルギーと対応づけるため)
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
        delta = - 2 * self.cell + 1
        return delta * self._get_site_de()

    def get_exchange_de(self, site_s0, site_s1):
        """
        spin 交換によるエネルギー差
        site_s0 は s = 0 のサイト
        site_s1 は s = 1 のサイトを指定
        site_s0, s1 の spin を交換した場合のエネルギー差を return する
        site_s0, s1 の表記は [x, y, z, site_id] の 4 要素の array
        最後に vertex をかける意味は
        spinが変わる点は vertex の回数だけ参照されるため
        othrer variables:
            tmp_cell: spin_flip 後の cell
        """
        tmp_cell = self.cell.copy()
        tmp_cell[tuple(site_s0.T)] = 1
        tmp_cell[tuple(site_s1.T)] = 0

        exchange_de = 0
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
            exchange_de += (dxi * eci) * vertex
        return exchange_de

    def get_exchange_de_small(self, site_s0, site_s1):
        """
        小さな cell 用の計算 (cell のサイズが相関関数のサイズより小さい場合)
        周期境界条件で eci がフリップするサイト自身を参照してしまうと
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

    def load_cell(self, fname):
        """
        size は揃えておく必要がある
        """
        with open(fname, 'rb') as rbfile:
            self.cell = pickle.load(rbfile)

    def get_conc(self):
        return self.cell.mean()


class CEMParser(object):
    """
    CEM の結果を読んで cluster, eci の情報を取り扱う

    Bcci に関しては対称操作が i-site 位置によって異なることが問題になる
    そこで clusters & ecis の組み合わせは各サイト毎の data として取り扱う
    したがって、 8 つの clusters, ecis のセットを作成する
    時間ができたら、他の class もこの方式に統一したい
    """

    @classmethod
    def get_unique_conv_site2idx(cls, clus, site_class):
        """
        対象操作を行って
        重複を削除 get_unique
        さらに、サイトの表記を site_class に基づいた
        index 表記に変換 (conv_site2idx)
        """
        unique_clus = []
        for i in range(len(clus[0])):
            sym_dups = cls._symm_cubic(clus[0][i], clus[1][i])
            sym_uni = cls._rm_dup(sym_dups)
            unique_clus.append(
                np.array([[site_class.conv_site2idx(site) for site in cluster]
                          for cluster in sym_uni]))
        return unique_clus

    @staticmethod
    def _symm_cubic(sites, equiv_pos=None):
        """
        対称操作で等価なサイトを作成する
        重複は考慮しない
        立方晶の method
        n体中の等価なサイトmにつき、どこを原点にするかという並進対称操作も必要
        m × 48 パターンを作成する
        副格子系は equiv_pos が必要
        args:
            sites: n体クラスターの座標 [n*3 of array]
            equiv_pos: 原点と等価なサイト位置 [m of list]
        """
        nbody = sites.shape[0]
        sites_t = []
        double = np.r_[sites, sites]
        # 先頭のサイトは常に[0, 0, 0]
        if not equiv_pos:
            equiv_pos = range(nbody)
        for i in equiv_pos:
            sites_t.append((double - sites[i])[i:i+nbody])
        sites_t = np.array(sites_t)

        symm_r = np.array([[0, 1, 2], [1, 0, 2], [0, 2, 1],
                           [1, 2, 0], [2, 0, 1], [2, 1, 0]])
        sites_r = []
        for site in sites_t:
            for symm in symm_r:
                sites_r.append(site[:, symm])
        sites_r = np.array(sites_r)

        symm_m = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1],
                           [-1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])
        return (sites_r * symm_m.reshape(8, 1, 1, 3)).reshape(48*len(equiv_pos),
                                                              nbody, 3)

    @classmethod
    def symm_cubic_bcci(cls, sites, site_class):
        """
        対称操作で等価なサイトを作成する
        重複は考慮しない
        立方晶の method
        n体中の等価なサイト m につき、どこを原点にするかという並進対称操作も必要
        等価なサイトは site_class から読み取る
        m × 48 パターンを作成する
        args:
            sites: n体クラスターの座標 [n*3 of array]
            equiv_pos: 原点と等価なサイト位置 [m of list]
        """
        nbody = sites.shape[0]
        symm_r = np.array([[0, 1, 2], [1, 0, 2], [0, 2, 1],
                           [1, 2, 0], [2, 0, 1], [2, 1, 0]])
        sites_r = []
        for symm in symm_r:
            for site in sites:
                sites_r.append(site[symm])
        sites_r = np.array(sites_r)
        symm_m = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1],
                           [-1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])
        sites_rm = (sites_r * symm_m.reshape(8, 1, 1, 3)).reshape(48, nbody, 3)
        equiv_label = [[site_class.EQUIV[site_class.conv_site2idx(site)[3]]
                        for site in cluster] for cluster in sites_rm]
        clusters = {}
        for label in set(BcciOctaSite.EQUIV.values()):
            tmp = []
            for s, eq in zip(sites_rm, equiv_label):
                for i in range(len(eq)):
                    double = np.r_[s, s]
                    if eq[i] == label:
                        tmp.append((double - s[i])[i:i+nbody])
            if tmp:
                # 原点も除去
                tmp = [x[1:] for x in cls._rm_dup(np.array(tmp))]
            clusters.update({label: tmp})
        return clusters

    @staticmethod
    def _symm_tetra(sites, equiv_pos=None, uniaxis='a'):
        """
        対称操作で等価なサイトを作成する
        重複は考慮しない
        正方晶の method
        n体のクラスターにつき、どこを原点にするかという並進対称操作も必要
        n体 × 16 パターンを作成する
        args:
            sites: n体クラスターの座標 [n*3 of array]
            uniax: 正方晶の 4 回対称軸方向を指定する
                   'a' or 'b' or 'c'
        """
        nbody = sites.shape[0]
        sites_t = []
        double = np.r_[sites, sites]
        # 先頭のサイトは常に[0, 0, 0]
        if not equiv_pos:
            equiv_pos = range(nbody)
        for i in equiv_pos:
            sites_t.append((double - sites[i])[i:i+nbody])
        sites_t = np.array(sites_t)

        symm_r = {'a': np.array([[0, 1, 2], [0, 2, 1]]),
                  'b': np.array([[0, 1, 2], [2, 1, 0]]),
                  'c': np.array([[0, 1, 2], [1, 0, 2]])}[uniaxis]
        sites_r = []
        for site in sites_t:
            for symm in symm_r:
                sites_r.append(site[:, symm])
        sites_r = np.array(sites_r)

        symm_m = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1],
                           [-1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])
        return (sites_r * symm_m.reshape(8, 1, 1, 3)).reshape(
            16*len(equiv_pos), nbody, 3)

    @staticmethod
    def _rm_dup(items):
        """
        重複する要素を除外する

        以下手順
        1. 第一カラム ([0, 0, 0]) を除去して list 化する
        2. 各要素をソートする
        3. 重複を除外する
        4. 第一カラムに [0, 0, 0] を戻して return

        結局除外するので 4. は不要かもしれないが、なんとなく気持ち悪いので入れておく
        """
        # 1.
        list_i = []
        for item in items[:, 1:, :]:
            li = [list(x) for x in item]
            list_i.append(li)
        # 2.
        sorted_items = [sorted(y) for y in list_i]
        # 3.
        unique = []
        for item in sorted_items:
            if not item in unique:
                unique.append(item)
        # 4.
        out = [[[0, 0, 0]] + x for x in unique]

        return out

    @classmethod
    def from_dirc_std(cls, dirc, site_class):
        """
        dirctory を指定して
        cluster の形 (site) を log.txt から、
        ecis を ecis.txt から読む
        表記法が conventional のケース
        """
        clus = cls.parse_logtxt_std(os.path.join(dirc, "log.txt"))
        ecis = cls.parse_ecitxt(os.path.join(dirc, 'eci.txt'))
        clus = [clus[i] for i in sorted(ecis.keys())]
        ecis = [ecis[i] for i in sorted(ecis.keys())]
        clus_out = [np.array([[site_class.conv_site2idx(site) for site in cl]
                              for cl in cluster]) for cluster in clus]
        # print(clus_out)
        with open(os.path.join(dirc, "cluster.pickle"), 'wb') as wbfile:
            pickle.dump({'clusters': clus_out[1:], 'ecis': ecis[1:],
                         'eci_point': ecis[0]}, wbfile)

    @classmethod
    def from_dirc_prim_fcc(cls, dirc, site_class):
        """
        dirctory を指定して
        cluster の形 (site) を log.txt から、
        ecis を ecis.txt から読む
        表記法が primitive のケース
        """
        clus = cls.parse_logtxt_prim_fcc(os.path.join(dirc, "log.txt"))
        ecis = cls.parse_ecitxt(os.path.join(dirc, 'eci.txt'))
        clus = [clus[i] for i in sorted(ecis.keys())]
        ecis = [ecis[i] for i in sorted(ecis.keys())]
        clus_out = [np.array([[site_class.conv_site2idx(site) for site in cl]
                              for cl in cluster]) for cluster in clus]
        with open(os.path.join(dirc, "cluster.pickle"), 'wb') as wbfile:
            pickle.dump({'clusters': clus_out[1:], 'ecis': ecis[1:],
                         'eci_point': ecis[0]}, wbfile)

    @classmethod
    def from_dirc_2sub(cls, dirc, site_class):
        """
        dirctory を指定して
        2 副格子のモデルを読み込む
        """
        clus, pos = cls.parse_logtxt_2sub(os.path.join(dirc, "log.txt"))
        ecis = cls.parse_ecitxt(os.path.join(dirc, 'eci.txt'))
        clus_dict = cls.extract_used_cluster(clus, pos, ecis)
        a = cls.get_unique_conv_site2idx(clus_dict['a'], site_class)
        c = cls.get_unique_conv_site2idx(clus_dict['c'], site_class)
        eci_a = clus_dict['a'][2]
        eci_c = clus_dict['c'][2]
        # print(a)
        with open(os.path.join(dirc, "cluster.pickle"), 'wb') as wbfile:
            pickle.dump({'sub_clusters': c, 'int_clusters': a,
                         'sub_ecis': eci_c, 'int_ecis': eci_a}, wbfile)

    @classmethod
    def from_dirc_bcci(cls, dirc, site_class):
        """
        dirctory を指定して
        bcci の構造を読み込む
        原点の [0, 0, 0, 0] 座標は抜いてある
        other variables:
            ecis: {id: ecis}
            clus: [座標の組み合わせの array]
            clusters: equiv_label で label 付けして原点を除去した cluster
        """
        clus = CEMParser.parse_logtxt_bcci(os.path.join(dirc, "log.txt"))
        clusters = []
        for c in clus:
            clusters.append(CEMParser.symm_cubic_bcci(c, BcciOctaSite))
        ecis = CEMParser.parse_ecitxt(os.path.join(dirc, "eci.txt"))
        out = {}
        for label in set(BcciOctaSite.EQUIV.values()):
            tmp_clus = []
            tmp_ecis = []
            for i in sorted(ecis.keys()):
                if clusters[i][label]:
                    tmp_clus.append(
                        np.array([[BcciOctaSite.conv_site2idx(site)
                                   for site in cluster]
                                  for cluster in clusters[i][label]]))
                    tmp_ecis.append(ecis[i])
            out.update({label:[tmp_clus, tmp_ecis]})
        with open(os.path.join(dirc, "cluster.pickle"), 'wb') as wbfile:
            pickle.dump(out, wbfile)

    @staticmethod
    def extract_used_cluster(cluster, clpos, ecis):
        """
        eci_ids に表記のある cluster のみを抜き出す
        """
        c_id = [i for i in sorted(ecis.keys()) if i in clpos[0]]
        a_id = [i for i in sorted(ecis.keys()) if i in clpos[2]]
        try:
            c_cls = [cluster[i] for i in c_id]
        except IndexError:
            print(len(cluster))
            print(c_id)
        c_ord = [clpos[1][i] for i in c_id]
        c_ecis = [ecis[i] for i in c_id]
        a_cls = [cluster[i] for i in a_id]
        a_ord = [clpos[3][i] for i in a_id]
        a_ecis = [ecis[i] for i in a_id]
        return {'c':[c_cls, c_ord, c_ecis], 'a':[a_cls, a_ord, a_ecis]}

    @staticmethod
    def meta_2sub(lines):
        """
        clusterの形状が記された位置 init - midle と
        clusterの元素種が記された位置 midle - end を re モジュールをつかって
        return する
        """
        meta_init = re.compile(r"number of clusters in DISORDERED STATE: .*")
        init = -1
        while not meta_init.match(lines[init]):
            init -= 1
        midle = init + 1
        meta_midle = re.compile(r"CF:     clust# clust#   dcrtn# #var- degen  dcrtn.*")
        while not meta_midle.match(lines[midle]):
            midle += 1
        end = midle + 1
        meta_end = re.compile(r"\*ECLI\*.*")
        while not meta_end.match(lines[end]):
            end += 1
        return init, midle, end

    @classmethod
    def parse_logtxt_2sub(cls, fname):
        """
        cluster の site を log.txt から読み込む
        i-s 用
        return:
            clus: 全ての cluster の site が入った array の list
            in_c: clus 中でラベル C の元素を含む list の番号
            pos_c: clus 中でラベル C の座標 (array) の番号
            in_a: clus 中でラベル A の元素を含む list の番号
            pos_a: clus 中でラベル A の座標 (array) の番号
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()

        init, midle, end = cls.meta_2sub(lines)

        clus = []
        for i in range(init+1, midle-1):
            if lines[i].split()[0] == 'cluster:':
                clus.append(np.array([[float(y) for y in x.split()]
                                      for x in lines[i+1:i+4]]).T)

        tmp = [line.split()[6] for line in lines[midle+2: end-4]]
        # print(tmp[0], tmp[-1]) # 最初と最後の pos が一致するか check
        pos_c = [[i for i in range(len(x)) if x[i] == 'C'] for x in tmp]
        in_c = [i for i in range(len(tmp)) if pos_c[i]]
        pos_a = [[i for i in range(len(x)) if x[i] == 'A'] for x in tmp]
        in_a = [i for i in range(len(tmp)) if pos_a[i]]

        return clus, (in_c, pos_c, in_a, pos_a)

    @classmethod
    def parse_logtxt_bcci(cls, fname):
        """
        bcci の site を log.txt から読み込む
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()

        init, midle, _ = cls.meta_2sub(lines)

        clus = []
        for i in range(init+1, midle-1):
            if lines[i].split()[0] == 'cluster:':
                clus.append(np.array([[float(y) for y in x.split()]
                                      for x in lines[i+1:i+4]]).T)
        return clus

    @classmethod
    def compare_pos_vs_clus(cls, fname, site_class):
        """
        sublattice の表記で
        log.txtから与えられた元素種と site 位置から判別した元素種とが一致するか
        check 用の method
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        _, midle, end = cls.meta_2sub(lines)
        log = [line.split()[6] for line in lines[midle+2: end-4]]
        clus = cls.parse_logtxt_bcci(fname)
        l = site_class.LABEL
        conv = ["".join([l[site_class.conv_site2idx(site)[3]]
                         for site in cluster])
                for cluster in clus]
        judge = []
        for label1, label2 in zip(log, conv):
            judge.append(label1 == label2)
        return judge

    @classmethod
    def parse_logtxt_std(cls, path):
        """
        log.txt を parse して cluster 群を return する
        表記方法が conventional cell 用
        以下手順
        1. ファイルの読み込み
        2. meta を使って cluster の情報が記載されている位置を習得
        3. parse
        4. 対称操作で cluster を複製、重複しないもののみを抜き出す
           また原点の座標は削除して return する
        """
        # 1.
        with open(path, 'r') as rfile:
            lines = rfile.readlines()
        # 2.
        meta_init = re.compile(r"number of clusters in DISORDERED STATE: .*")
        init = 0
        while not meta_init.match(lines[init]):
            init += 1
        end = init + 1
        meta_end = re.compile(r"Species represented in concentration vector: .*")
        while not meta_end.match(lines[end]):
            end += 1
        # 3.
        tmp_clus = []
        for i in range(init+1, end-1):
            if lines[i].split()[0] == 'cluster:':
                tmp_clus.append(np.array([[float(y) for y in x.split()]
                                          for x in lines[i+1:i+4]]).T)
        # 4.
        clus = [[x[1:] for x in cls._rm_dup(cls._symm_cubic(sub_clus))]
                for sub_clus in tmp_clus]

        return clus

    @classmethod
    def parse_logtxt_prim_fcc(cls, path):
        """
        log.txt を parse して cluster 群を return する
        表記方法が primitive fcc cell 用
        以下手順
        1. ファイルの読み込み
        2. meta を使って cluster の情報が記載されている位置を習得
        3. parse
        4. デカルト座標に変換
        5. 対称操作で cluster を複製、重複しないもののみを抜き出す
           また原点の座標は削除して return する
        """
        # 1.
        with open(path, 'r') as rfile:
            lines = rfile.readlines()
        # 2.
        meta_init = re.compile(r"number of clusters in DISORDERED STATE: .*")
        init = 0
        while not meta_init.match(lines[init]):
            init += 1
        end = init + 1
        meta_end = re.compile(r"Species represented in concentration vector: .*")
        while not meta_end.match(lines[end]):
            end += 1
        # 3.
        prim = []
        for i in range(init+1, end-1):
            if lines[i].split()[0] == 'cluster:':
                prim.append(np.array([[float(y) for y in x.split()]
                                      for x in lines[i+1:i+4]]).T)
        axis = np.array([[1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2]])
        # 4.
        tmp_clus = [np.array([(axis * cood.T).sum(-1).T for cood in clus])
                    for clus in prim]
        # print(prim)
        # print(tmp_clus)

        # 5.
        clus = [[x[1:] for x in cls._rm_dup(cls._symm_cubic(sub_clus))]
                for sub_clus in tmp_clus]

        return clus

    @staticmethod
    def conv_NaCl(site):
        """
        NaCl site の表記に変換する
        dcimal to integer
        """
        integ = []
        decim = []
        for coord in site:
            integ.append(math.floor(coord))
            decim.append(coord - math.floor(coord))
        site_id = {(0.0, 0.0, 0.0): 0, (0.0, 0.5, 0.5): 1,
                   (0.5, 0.0, 0.5): 2, (0.5, 0.5, 0.0): 3,
                   (0.5, 0.0, 0.0): 4, (0.0, 0.5, 0.0): 5,
                   (0.0, 0.0, 0.5): 6, (0.5, 0.5, 0.5): 7}[tuple(decim)]
        integ.append(site_id)
        return integ

    @staticmethod
    def parse_ecitxt(fname):
        """
        null クラスターは flip に関与しないので削除する
        clus_pos と対応させるために id は -1 する
        return:
            {id: eci}
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        # num_cls = int(lines[0].split()[0])
        pos = lines.index(" 0.000000\n")
        cls_ids = []
        for i in range(1, pos):
            cls_ids += [int(x) - 1 for x in lines[i].split()]
        ecis = []
        for i in range(pos+1, len(lines)):
            ecis += [float(x) for x in lines[i].split()]
        ecis_dict = {i:x for i, x in zip(cls_ids, ecis)}
        ecis_dict.pop(-1)
        return ecis_dict


class CEMParser02(object):
    """
    CEM を parse するための method
    bcci の手法を基礎にして整理し直す
    """
    @staticmethod
    def parse_ecitxt(fname):
        """
        null クラスターは flip に関与しないので削除する
        clus_pos と対応させるために id は -1 する
        return:
            {id: eci}
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        # num_cls = int(lines[0].split()[0])
        pos = lines.index(" 0.000000\n")
        cls_ids = []
        for i in range(1, pos):
            cls_ids += [int(x) - 1 for x in lines[i].split()]
        ecis = []
        for i in range(pos+1, len(lines)):
            ecis += [float(x) for x in lines[i].split()]
        ecis_dict = {i:x for i, x in zip(cls_ids, ecis)}
        ecis_dict.pop(-1)
        return ecis_dict

    @classmethod
    def symm_cubic(cls, sites, site_class):
        """
        対称操作で等価なサイトを作成する
        重複は考慮しない
        立方晶の method
        n体中の等価なサイト m につき、どこを原点にするかという並進対称操作も必要
        等価なサイトは site_class から読み取る
        m × 48 パターンを作成する
        args:
            sites: n体クラスターの座標 [n*3 of array]
            equiv_pos: 原点と等価なサイト位置 [m of list]
        """
        nbody = sites.shape[0]
        symm_r = np.array([[0, 1, 2], [1, 0, 2], [0, 2, 1],
                           [1, 2, 0], [2, 0, 1], [2, 1, 0]])
        sites_r = []
        for symm in symm_r:
            for site in sites:
                sites_r.append(site[symm])
        sites_r = np.array(sites_r)
        symm_m = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1],
                           [-1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])
        sites_rm = (sites_r * symm_m.reshape(8, 1, 1, 3)).reshape(48, nbody, 3)
        equiv_label = [[site_class.EQUIV[site_class.conv_site2idx(site)[3]]
                        for site in cluster] for cluster in sites_rm]
        clusters = {}
        for label in set(site_class.EQUIV.values()):
            tmp = []
            for s, eq in zip(sites_rm, equiv_label):
                for i in range(len(eq)):
                    double = np.r_[s, s]
                    if eq[i] == label:
                        tmp.append((double - s[i])[i:i+nbody])
            if tmp:
                # 原点も除去
                tmp = [x[1:] for x in cls.rm_dup(np.array(tmp))]
            clusters.update({label: tmp})
        return clusters

    @staticmethod
    def rm_dup(items):
        """
        重複する要素を除外する

        以下手順
        1. 第一カラム ([0, 0, 0]) を除去して list 化する
        2. 各要素をソートする
        3. 重複を除外する
        4. 第一カラムに [0, 0, 0] を戻して return

        結局除外するので 4. は不要かもしれないが、なんとなく気持ち悪いので入れておく
        """
        # 1.
        list_i = []
        for item in items[:, 1:, :]:
            li = [list(x) for x in item]
            list_i.append(li)
        # 2.
        sorted_items = [sorted(y) for y in list_i]
        # 3.
        unique = []
        for item in sorted_items:
            if not item in unique:
                unique.append(item)
        # 4.
        out = [[[0, 0, 0]] + x for x in unique]

        return out

    @classmethod
    def from_dirc_bcci(cls, dirc, site_class):
        """
        dirctory を指定して
        bcci のクラスターの情報を読み込む
        原点の [0, 0, 0, 0] 座標は抜いてある
        other variables:
            ecis: {id: ecis}
            clus: [座標の組み合わせの array]
            clusters: equiv_label で label 付けして原点を除去した cluster
        """
        clus = cls.parse_logtxt(os.path.join(dirc, "log.txt"))
        clusters = []
        for c in clus:
            clusters.append(cls.symm_cubic(c, site_class))
        ecis = cls.parse_ecitxt(os.path.join(dirc, "eci.txt"))
        out = {}
        for label in set(site_class.EQUIV.values()):
            tmp_clus = []
            tmp_ecis = []
            for i in sorted(ecis.keys()):
                if clusters[i][label]:
                    tmp_clus.append(
                        np.array([[site_class.conv_site2idx(site)
                                   for site in cluster]
                                  for cluster in clusters[i][label]]))
                    tmp_ecis.append(ecis[i])
            out.update({label:[tmp_clus, tmp_ecis]})
        with open(os.path.join(dirc, "cluster.pickle"), 'wb') as wbfile:
            pickle.dump(out, wbfile)

    @classmethod
    def parse_logtxt(cls, fname):
        """
        site を log.txt から読み込む
        ToDo: _meta_2sub が 2 sublattice 以外も使えるか check する
              その上で midle と end の両方は不要かもしれない
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()

        init, midle, _ = cls._meta_2sub(lines)

        clus = []
        for i in range(init+1, midle-1):
            if lines[i].split()[0] == 'cluster:':
                clus.append(np.array([[float(y) for y in x.split()]
                                      for x in lines[i+1:i+4]]).T)
        return clus

    @staticmethod
    def _meta_2sub(lines):
        """
        clusterの形状が記された位置 init - midle と
        clusterの元素種が記された位置 midle - end を re モジュールをつかって
        return する
        """
        meta_init = re.compile(r"number of clusters in DISORDERED STATE: .*")
        init = -1
        while not meta_init.match(lines[init]):
            init -= 1
        midle = init + 1
        meta_midle = re.compile(r"CF:     clust# clust#   dcrtn# #var- degen  dcrtn.*")
        while not meta_midle.match(lines[midle]):
            midle += 1
        end = midle + 1
        meta_end = re.compile(r"\*ECLI\*.*")
        while not meta_end.match(lines[end]):
            end += 1
        return init, midle, end
