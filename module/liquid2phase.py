#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
液体の自由エネルギーを"2 phase model"を使って算出する
[1] Zi-Kui Liu et al. Fluid Phase Equilibria, 360, 44 (2013)
[2] Norman F. Carnahan & Kenneth E. Starling, J. Chem. Phys., 53, 600 (1970)

class
    VDoS: vibrational density of state
          freaquancy [THz]
          DoS [1/THz] のデータを格納
"""
from __future__ import division
import os
import pylab
import pickle
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

__date__ = "Jan 28 2015"


def main():
    """main"""
    #Al30Fe()
    # Al()
    # Fe()
    return
    TEST_PATH = ("/Users/enoki/Researches/Analysis/Codes/01_testRun/vdos")
    # f = os.path.join(PATH, 'VDOS.txt')
    # a = VDoS.from_file(f, 2300, 200, 4311.57, 26.9815395, header=1)
    f = os.path.join(TEST_PATH, 'AlCu256/vdos.pickle')
    a = VDoS.from_pickle(f, 2300, 256*2, 4311.57, 26.9815395)
    print(a._integral(a.vdos))
    print(a.get_entropy_gas())
    print(a.get_entropy_solid())
    gas = a.get_vdos_gas()(a.vdos[:, 0])
    pylab.plot(a.vdos[:, 0], a.vdos[:, 1])
    pylab.plot(a.vdos[:, 0], gas)
    pylab.plot(a.vdos[:, 0], a.get_vdos_solid())
    pylab.show()
    # b = VDoS.from_file(f, 1900, 200, 4208.45, 26.9815395, header=1)
    # b._s0 = 82e-12
    # print(b.get_entropy_gas())

def Al():
    """
    Al
    srcファイル(Al.txt)はテストディレクトリに移した
    01_testRun/vdos/liquid2phase
    ToDo:そのうちtestモジュールを作成する
    """
    fname = ("/Users/enoki/Desktop/Al.txt")
    with open(fname, 'r') as rfile:
        read_lines = rfile.readlines()
    Al = np.array([[float(x) for x in line.split()[0:2]]
                     for line in read_lines[1:]])[:, 0:2]

    pylab.plot(Al[:,0], Al[:,1])

    al = VDoS(Al, 2300, 200, 4311.57, 26.9815395)
    print(al.get_entropy_gas())
    print(al.get_entropy_solid())
    print(al._Delta)
    print(al._f)

    print(al._integral(Al))

    pylab.show()

def Fe():
    """
    Fe 鉄がなぜ合わないのか
    横軸のスケールを誤っている可能性がある
    """
    fname = ("/Users/enoki/Desktop/Fe.txt")
    with open(fname, 'r') as rfile:
        read_lines = rfile.readlines()
    Fe = np.array([[float(x) for x in line.split()[0:2]]
                     for line in read_lines[1:]])[:, 0:2]

    pylab.plot(Fe[:,0], Fe[:,1])
    Fe[:,0] *= (1/2*3)
    # Fe[:,1] *= 300/326

    fe = VDoS(Fe, 1900, 100, 1147.69, 55.847)
    print(fe.get_entropy_gas())
    print(fe.get_entropy_solid())
    print(fe._Delta)
    print(fe._f)

    print(fe._integral(Fe))

    gas = fe.get_vdos_gas()(fe.vdos[:, 0])
    pylab.plot(fe.vdos[:, 0], gas)
    pylab.plot(fe.vdos[:, 0], fe.get_vdos_solid())

    pylab.show()

def Al30Fe():
    """
    Al30Fe
    VDoSは積分値は 3N
    各元素毎に分けて、それぞれを規格化する 3Na + 3Nb +...
    この規格化したVDoSから、それぞれの Δ, f が求まる
    エントロピーに mass の fraction をかけて合算する

    文献はこの手順で計算しているのだが、
    本当にそれでいいのだろうか？？
    """
    fname = ("/Users/enoki/Desktop/dos_Al-30Fe.txt")
    with open(fname, 'r') as rfile:
        read_lines = rfile.readlines()
    Al = np.array([[float(x) for x in line.split()[0:2]]
                     for line in read_lines[1:]])[:, 0:2]

    Fe = np.array([[float(x) for x in line.split()[2:4]]
                     for line in read_lines[1:22]])[:, 0:2]
    Fe[:, 1] *= 25 / 53.6306
    Fe[:, 0] = 40 - Fe[:, 0]

    pylab.plot(Al[:,0], Al[:,1])
    pylab.plot(Fe[:,0], Fe[:,1])

    ave = (55.847 * 0.6 + 26.9815395 * 1.4) / 2

    al = VDoS(Al, 1900, 140, 3089.81, 26.9815395)
    fe = VDoS(Fe, 1900, 60, 3089.81, 55.847)

    print('Al')
    print(al.get_entropy_gas())
    print(al.get_entropy_solid())
    print(al._Delta)
    print(al._f)

    print('Fe')
    print(fe.get_entropy_gas())
    print(fe.get_entropy_solid())
    print(fe._Delta)
    print(fe._f)

    print('total')
    print((al.get_entropy_gas()*0.7*26.9 + fe.get_entropy_gas()*0.3*55.847)/ave)
    print((al.get_entropy_solid()*0.7*26.9 + fe.get_entropy_solid()*0.3*55.847)/ave)

    print(al._integral(Al)[0])
    print(al._integral(Fe)[0])
    f1 = interp1d(Al[:, 0], Al[:, 1], kind='cubic')
    f2 = interp1d(Fe[:, 0], Fe[:, 1], kind='cubic')
    x = np.arange(0.1, 39, 0.1)
    y = [f1(i) + f2(i) for i in x]

    pylab.plot(x, y)





    pylab.show()


def get_s_con(n, x, f):
    """
    configuration の etropy を retun
    """
    kb = 8.31446
    sigma_nixi = (n * x).sum()
    nixifij = (n * x).reshape(-1, 1) * f
    print(nixifij.sum(axis=0).reshape(1, -1) / nixifij)
    ln = np.log(nixifij.sum(axis=0).reshape(1, -1) / nixifij)
    sigma = nixifij * ln

    print(nixifij)
    print(ln)
    print(sigma.sum() / sigma_nixi * kb)

# n = np.array([11*2, 9.5*2])
# x = np.array([0.8, 0.2])
# f = np.array([[6/11, 5/11], [0.99, 0.01]])

# get_s_con(n, x, f)






class VDoS(object):
    """
    vibrational density of state
    attributes
        vdos: np.array n × 2
              vdos[:, 0]: freaquancy [THz]
              vdos[:, 1]: DoS [1/THz]
        temp: temperature [K]
        num_atoms: number of atoms in cell
        volume: volume of cell [Å^3]
        mass: atomic weight (g/mol)
        d0: VDoS at v=0 (1/THz)
    """
    def __init__(self, vdos, temp, num_atoms, volume, mass):
        self._kb = 1.3806488e-23
        self._nA = 6.023e23
        self._planck = 6.62606957e-34 # planck constant [Js]

        self.vdos = vdos
        self._m = mass / self._nA / 1000 # [kg]
        self._beta = 1 / self._kb / temp # [1/J]
        self._d0 = vdos[0, 1] * 1e-12 # [1/s]
        self._v = volume * 1e-30 # [m^3]
        self._n = num_atoms
        self._Delta = self.get_dimless_diffus()
        self._f, _ = self.get_liq_frac(self._Delta)
        self._y = (self._f / self._Delta) ** (2. / 3)
        self.vdos_sol = self.get_vdos_solid()
        self.vdos_gas = self.get_vdos_gas()(self.vdos[:, 0])
        self.w_gas = self.get_weight_func_gas()
        # function
        self.w_sol = self.get_weight_func_solid()

    @classmethod
    def from_pickle(cls, fname, temp, num_atoms, volume, mass):
        """
        vdos の pickle から load する
        data の format は
        1 列目: freaquancy [THz]
        2 列目: DoS [THz]
        """
        with open(fname, 'rb') as rbfile:
            data = pickle.load(rbfile)
        return cls(data, temp, num_atoms, volume, mass)

    @classmethod
    def from_file(
        cls, fname, temp, num_atoms, volume, mass, scale=1, header=0):
        """
        ファイルからVDoSのデータ読み込み
        ファイルのフォーマットは
        1 列目: freaquancy [THz]
        2 列目: DoS [THz]
        またheaderの行数分スキップしてデータを読む
        単位の換算が必要な場合 scaleに換算値を入れる
        freaquancy * scale
        DoS / scaleが読み込まれる
        """
        with open(fname, 'r') as rfile:
            read_lines = rfile.readlines()
        data = np.array([[float(x)/scale for x in line.split()]
                         for line in read_lines[header:]])[:, 0:2]
        return cls(data, temp, num_atoms, volume, mass)

    @staticmethod
    def _integral(data, kind='cubic'):
        """
        (x, y)の離散データを線形補間して積分値を求める
        補間の関数形はkindで指定
        積分の範囲はx_minからx_maxまで
        dataはxの昇順にソートしておく
        """
        f = interp1d(data[:, 0], data[:, 1], kind=kind)
        return quad(f, data[0, 0], data[-1, 0])

    @staticmethod
    def equation_3(delta, frac):
        """
        ref[1] eq(3) の Δ と f の関係式における左
        """
        return (2 * delta ** (-9./2) * frac ** (15./2) -
                6 * delta ** (-3) * frac ** 5 -
                delta ** (-3. / 2) * frac ** (7. / 2) +
                6 * delta ** (-3. / 2) * frac ** (5. / 2) +
                2 * frac - 2)

    @classmethod
    def get_liq_frac(cls, delta, initial_frac=0, accuracy=1e-6):
        """
        liquid fraction f
        解析的に解けないので(maybe)eq(3)を使って摺り合わせる
        """
        frac = initial_frac
        residual = cls.equation_3(delta, frac)
        while residual < 0:
            frac += accuracy
            residual = cls.equation_3(delta, frac)
        return frac, residual


    def get_dimless_diffus(self):
        """
        dimension less diffusion constant Δ
        ref[1] eq(4)
        """
        return (2 * self._d0 / 9 / self._n *
                (np.pi / self._beta/ self._m) ** (0.5) *
                (self._n / self._v) ** (1. / 3) *
                (6 / np.pi) ** (2. / 3))

    def get_vdos_gas(self):
        """
        function of vdos of gas compornent
        unit is [1/THz]
        ref[1] eq(2)
        """
        def vdos(v):
            return (
                self._d0 / (1 + ((np.pi * self._d0 * v * 1e12) /
                            (6 * self._f * self._n)) ** 2 )) * 1e12
        return vdos

    def get_vdos_solid(self):
        """
        vdos of solid compornent
        unit is [1/THz]
        ref[1] eq(1)
        """
        vdos_gas = self.get_vdos_gas()
        return self.vdos[:, 1] - vdos_gas(self.vdos[:, 0])

    def get_weight_func_solid(self):
        """
        weight function of solid w_sol(v)
        ref[1] eq(6)
        """

        def w_sol(v):
            bhv = self._beta * self._planck * v * 1e12  # v [1/THz] -> [1/Hz]
            bhv += 1e-15
            return bhv / (np.exp(bhv) - 1) - np.log(1 - np.exp(-bhv))
        return w_sol

    def get_entropy_solid(self):
        """
        vibrational entropy of solid compornent
        ref[1] eq(5) second term
        """
        weighted_vdos = (self.vdos_sol * self.w_sol(self.vdos[:, 0]) *
                         self._kb * self._nA / self._n)
        data = np.c_[self.vdos[:, 0], weighted_vdos]
        return self._integral(data)[0]

    def get_weight_func_gas(self):
        """
        weight function of gas compornent
        ref[1] eq(7)

        理想気体のエントロピーがvに依存しない場合
        Tと圧力のみに依存する (vに依らずconstant になる)
        v に依存する形で表現できるのかもしれない
        しかし、vdos_gas の v 依存性は形が決まっているので、
        結局、積をとって積分すると用いている値になるのではないかと考えられる (未確認)
        """
        c_t = (2 * np.pi * self._m / self._beta / self._planck ** 2) ** 1.5
        # entropy of ideal gas [per atom]  devided by kb
        s_ideal_gas = (np.log(self._v * c_t / self._n) + 2.5)
        fy = self._f * self._y
        return 1. / 3 * (s_ideal_gas +
                         np.log(1 + fy + fy ** 2 - fy ** 3) / (1 - fy) ** 3 +
                         (fy * (3 * fy - 4)) / (1 - fy) ** 2)

    def get_entropy_gas(self):
        """
        vibrational entropy of gas compornent
        [1] eq(5) first term
        積を積分するとentropyが小さく評価される
        原因は積分の収束が遅く積分を十分高い v まで実行する必要があるためとわかった

        weightがconstantだとすると weight は積分の外に出せる
        また vdos_gas の積分値は f×3×N に収束する
        積分だと収束が遅いのでこの値を weight に掛けて算出する
        """
        return self._f * 3 * self.w_gas * self._kb * self._nA

if __name__ == "__main__":
    main()




