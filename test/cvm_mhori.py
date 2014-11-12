#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
自然逐次計算法によるcvmコード
(材料システム工学)
e_ijとeci(=v_i)は違う
このコードではクラスターのエネルギーeijを使う
相関関数は使用せずクラスター濃度x_i, y_ijを使用
"""
import itertools
import copy
import math

class Configration(object):
    """
    配列を作成
    originalのM行列に対応
    """
    @classmethod
    def mat_clus_tetra(self):
        config = [format(int(self.base10to(x,2)), '0>4') for x in range(16)]
        def point(i):
            return int(i, 2)
        def pair(i):
            return int(i, 2) + 2
        def get_m(i):
            m = []
            m += [point(x) for x in i]
            combi = list(itertools.combinations('0123', 2))
            m += [pair(i[int(x[0])]+i[int(x[1])]) for x in combi]
            return m
        mat_clus = [get_m(x) for x in config]
        return mat_clus

    @classmethod
    def base10to(cls, n, b):
        if (int(n/b)):
            return cls.base10to(int(n/b), b) + str(n%b)
        return str(n%b)


NIT = Configration.mat_clus_tetra()  # 16種のconfiguration matrix
ETA = 1.e-9  # 収束の精度
NTC = 3  # size of ttk and ccp
TTK = [1.0, 1.1, 1.1]  # temperatures
CCP = [0.0, 0.0, 0.5]  # effective chemical potentials (A原子とB原子の差分)
EIJ = [0.0] * 30  # 30種のクラスターのenergy
EIJ[2:5] = [1.0, -1.0, -1.0, 1.0]  # 2体以外のクラスターのエネルギーは0とする
XI = [0.5] * 30  # point-tetraまでの30種のクラスターの出現確率(x_i, y_ij,...)

def nat_itr(cp, eij, eta, nit, tk, xi, init):
    ms = -1  # minus spin
    xinit = copy.deepcopy(xi)  # xiの初期値を保管 convergenceの判定に使う

    n = 13
    zch = [0] * 30  # 計算用の容器 (30種のクラスターに対応)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    i7 = ms ** i
                    j7 = ms ** j
                    k7 = ms ** k
                    l7 = ms ** l
                    n += 1
                    nnn = n - 14
                    sxy = [xi[x] for x in nit[nnn]]  # 2体までの出現確率
                    pp = float(i7 + j7 + k7 + l7)
                    p = pp * cp / 4.0

                    # クラスター変分法 (10.28)式1
                    # 四面体クラスターのエネルギー　=> 二体クラスターを数えてsumをとる
                    eijkl = sum([eij[x] for x in nit[nnn][4:10]])
                    # chemical potentialを引いたエネルギー
                    # -e4/2βが(10.28)のexpの中身
                    e4 = eijkl - p

                    # 積をとる (logを利用して計算)
                    x12 = sum([math.log(sxy[x]) for x in range(4)])
                    dex12 = math.exp(x12)

                    y22 = sum([math.log(sxy[x]) for x in range(4,10)])
                    dey22 = math.exp(y22)

                    # クラスター変分法 (4.23)式
                    zch[n] = (math.exp(e4/(-2*tk))*(dey22**0.5) /
                              (dex12**(5.0/8.0)))

    ztt = sum(zch[14:30])
    # ラグランジュ乗数 => 格子点あたりのグランドポテンシャル
    gtmp = -2.0*tk*math.log(ztt)
    print(gtmp)
    for i in range(14, 30):
        zch[i] *= math.exp(gtmp/(2.0*tk))

    # クラスター変分法 (4.24)式 (縮小の式)
    xi[0:6] = [0.0]*6
    xi[0] += sum(zch[14:22])
    xi[1] += sum(zch[22:30])
    xi[2] += sum(zch[14:18])
    xi[3] += sum(zch[18:22])
    xi[4] += sum(zch[22:26])
    xi[5] += sum(zch[26:30])
    xi[6] = zch[14] + zch[15]
    xi[7] = zch[16] + zch[17]
    xi[8] = zch[18] + zch[19]
    xi[9] = zch[20] + zch[21]
    xi[10] = zch[22] + zch[23]
    xi[11] = zch[24] + zch[25]
    xi[12] = zch[26] + zch[27]
    xi[13] = zch[28] + zch[29]
    xi[14:30] = zch[14:30]

    initf = init
    # A value to examine the convergency　誤差
    dxnit = sum([math.fabs(xinit[x]-xi[x]) for x in range(14)])
    return dxnit




i = 2
for init in range(2000):
    dxnit = nat_itr(CCP[i], EIJ, ETA, NIT, TTK[i], XI, init)
    if dxnit < ETA:
        break

print(XI)
