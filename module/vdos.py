#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
VDOSを計算する
"""
import os
import numpy as np
import errno
import pickle
from scipy import signal
from scipy.integrate import quad
from scipy.interpolate import interp1d
import pylab

__date__ = "Apr 3 2015"

def main():
    """
    main
    """
    pass


class Atoms(object):
    """
    各時間での原子の座標データを保管
    速度 vacf, vdos を計算
    XDATCAR を読んで生成
    """
    def __init__(self, elements, num_atoms, coods, dt, head):
        self.elements = elements
        self.num_atoms = num_atoms
        self.coods = coods
        self.velocity = self.get_velocity_from_coods(coods)
        self.dt = dt
        self.POMAS = {"Cu": 63.546, "Al": 26.982, "Fe": 55.847}
        self.head = head

    def make_poscars(self, num, dst="."):
        """
        MD 計算における enthalpy 計算用の method
        coods を遡って POSCAR を num の数だけ作成する
        """
        for i in range(num):
            dir_name = "pos_{0:03d}".format(i)
            path = os.path.join(dst, "enthalpy", dir_name)
            os.makedirs(path, exist_ok=True)
            lines = "".join(self.head)
            cood_list = ["  ".join([str(x) for x in y]) for y in self.coods[-i]]
            tmp = "\n".join(cood_list)
            lines += "Direct\n"
            lines += tmp
            path = os.path.join(path, "POSCAR")
            with open(path, 'w') as wfile:
                wfile.write(lines)

    def get_volume(self):
        """
        head 読んだから volume を return
        対角行列のみ対応
        """
        scale = float(self.head[1])
        mat = np.array([[float(x) for x in line.split()] for line in self.head[2:5]])
        return scale * mat[0, 0] * mat[1, 1] * mat[2, 2]

    @staticmethod
    def _mkdir(path):
        """
        Alternative command of 'mkdir -p'
        At Python3.3, we enable to use os.makedirs(path, exist_ok=True)
        """
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise


    def get_cn_single(self, i, j, rmin, pos, d_mesh):
        """
        coodination number を return する
        グラフから式 22 の Rmin を決定する
        1. 各座標の cmatrix 作成 (原子数, 原子数, 座標)
        2. 差分 diff を取る
        3. 0.5 を閾値に mod を取る
        4. 距離 dis を計算
        5. rmin を区切った mesh を作成
        6. dis を mesh 分に複製
        7. nij を計算 (式22)
        8. dnij を計算
        """
        # 1.
        coods = self.coods[pos]
        istart = sum(self.num_atoms[:i])
        iend = istart + self.num_atoms[i]
        jstart = sum(self.num_atoms[:j])
        jend = jstart + self.num_atoms[j]
        coods_j = coods[jstart:jend]
        tmp = np.arange(istart, iend).reshape(-1, 1)
        zero = np.zeros_like(np.arange(len(coods_j))).reshape(1, -1)
        idx = tmp + zero
        cmatrix = coods[idx]
        # 2.
        diff = coods_j - cmatrix
        # 3.
        diff[diff > 0.5] = 1 - diff[diff > 0.5]
        diff[diff < -0.5] = 1 + diff[diff < -0.5]
        # 4.
        dis = np.sqrt((diff * diff).sum(-1))
        # 5.
        mesh = np.arange(0, 1, d_mesh).reshape(-1, 1, 1)
        # 6.
        zero = np.zeros_like(np.arange(len(mesh))).reshape(-1, 1, 1)
        dis_mesh = zero + dis
        # 7.
        nij = (dis_mesh < mesh).sum(-1).mean(-1)
        # 8.
        dnij = nij[1:] - nij[:-1]
        ddnij = dnij[1:] - dnij[:-1]
        dr = np.arange(0, (len(ddnij)+1)*d_mesh, d_mesh)
        return nij, dnij, dr
        # pylab.plot(nij)
        pylab.plot(dr, dnij)
        pylab.plot([rmin, rmin], [0, 5])
        pylab.show()
        print(nij)

    def get_cn(self, rmin, d_mesh):
        cutoff = int(rmin//d_mesh)+1
        num = len(self.num_atoms)
        frac = np.zeros_like(np.arange(0, num**2, 1.)).reshape(num, num)
        n_out = []
        for i in range(2):
            ni = 0
            for j in range(2):
                nij, dnij = 0, 0
                for pos in range(100):
                    n, d, r = self.get_cn_single(i, j, rmin, pos, d_mesh)
                    nij += n
                    dnij += d
                    ni += n
                print(nij[cutoff] / 100.)
                frac[i][j] = nij[cutoff] / 100.
                print(frac)
                # print(i, j)
                # print(nij[cutoff] / 100)
                label = "{0} {1}".format(i, j)
                pylab.plot(r, dnij / 100, label=label)
            frac[i, :] /= ni[cutoff] / 100
            # print("ni", i)
            n_out.append(ni[cutoff] / 100)
            # print(ni[cutoff] / 100)
        pylab.legend(loc="upper left")
        print(r[cutoff])
        # pylab.plot([rmin, rmin], [0, 5])
        # pylab.savefig('/Users/enoki/Desktop/plot.eps')
        pylab.show()
        n_out = np.array(n_out)
        return(n_out, frac)

    @staticmethod
    def get_velocity_from_coods(coods):
        """
        velocity を return する
        """
        velocity = coods[1:] - coods[0: -1]
        velocity[velocity > 0.9] = 1 - velocity[velocity > 0.9]
        velocity[velocity < -0.9] = 1 + velocity[velocity < -0.9]
        return velocity

    @classmethod
    def from_xdatcar(cls, fname='XDATCAR', dt=5.0, start=0, end=-1):
        """
        XDATCAR を読む
        dt は INCAR で設定した POTIM を入力
        単位は fs
        周期境界条件で座標が飛ぶ点があるため、 velocity は補正する
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        elements = lines[5].split()
        num_atoms = [int(x) for x in lines[6].split()]
        total_atoms = sum(num_atoms)
        steps = range(int((len(lines) - 7) / (total_atoms + 1)))
        coods = np.array(
            [[[float(x) for x in line.split()]
              for line in lines[8+i*(total_atoms+1):
                                8+i*(total_atoms+1)+total_atoms]]
             for i in range(steps[start], steps[end])])
        head = lines[0:7]
        return cls(elements, num_atoms, coods, dt, head)

    def get_vacf(self):
        """
        速度自己相関関数を return する
        half: データの半分を使う
        a: (1, half, 原子数, 3vector)
        vacf_zero: dt=0 の速度自己相関関数 (規格化に用いる)
        decay: dt を取り扱う matrix (half, half)
        vacf: (dt, 原子数)
        """
        half = len(self.velocity)//2
        a = self.velocity[:half, :, :].reshape(1, half, -1, 3)
        arb_a = (a * a).sum(axis=-1).sum(axis=1)
        decay = np.arange(half).reshape(-1, 1) + np.arange(half).reshape(1, -1)
        b = self.velocity[decay, :, :]
        arb_b = (b * b).sum(axis=-1).sum(axis=1)
        vacf = (a * b).sum(axis=-1).sum(axis=1) / np.sqrt(arb_a * arb_b)
        return vacf

    def get_vdos(self):
        """
        1. 終端を 0 に補正する dump
        2. 横軸を dt から求める [THz]
        3. vdos の積分値は 3 * num_atoms に規格化
           (vacf の横軸の L は len(vacf) * dt)
           (それぞれの成分で規格化を行う)
        4. element 毎でラベル付けした dict で return する
        """
        vacf = self.get_vacf()
        mirror_dat = vacf[::-1]
        vacf_symm = np.r_[vacf, mirror_dat]
        # window = 1 - np.hanning(len(vacf_symm)).reshape(-1, 1)
        # window の幅で D0 などの値が結構変わる...
        # window の取り方に任意性が残る
        # アインシュタインの方法と D を比較することは重要かもしれない
        # ただ幸いなことに、 S_v (Sgas + Ssol) は
        # あまり window に依存していないように見える
        # window = signal.gaussian(len(vacf_symm), std=len(vacf_symm)/50).reshape(-1, 1)
        window = signal.gaussian(len(vacf_symm), 20).reshape(-1, 1)
        window = np.r_[window, window][len(vacf_symm)/2:1.5*len(vacf_symm)]
        vacf_symm *= window
        i = 0
        data = {}
        for num, elem in zip(self.num_atoms, self.elements):
            # waited_vacf = vacf_symm[:, i:i+num].sum(axis=-1) * self.POMAS[elem]
            # vdos = np.fft.fft(waited_vacf)[:vacf.shape[0]]
            vdos = np.fft.fft(vacf_symm[:, i:i+num].sum(axis=-1))[:vacf.shape[0]]
            xaxis = np.arange(0, vdos.shape[0], 1.0) / \
                 self.dt / 1e-3 / vdos.shape[0] / 2
            data.update({elem: vdos.real})
            i += num
            integ = self._integral(xaxis, vdos.real)
            col = 3 * num / integ[0]
            data[elem] *= col

        return xaxis, data

    def get_averaged_speed(self):
        """
        平均の速さをreturn
        """
        return ((self.velocity ** 2) ** 0.5).mean()

    @staticmethod
    def _integral(datax, datay, kind='cubic'):
        """
        datax, datay の離散データを線形補間して積分値を求める
        補間の関数形はkindで指定
        積分の範囲はx_minからx_maxまで
        dataはxの昇順にソートしておく
        """
        f = interp1d(datax, datay, kind=kind)
        return quad(f, datax[0], datax[-1])

    def save_pickle(self, fname):
        """
        vdos を pickle に保存
        """
        xaxis, vdos = self.get_vdos()
        for elem in self.elements:
            data = np.c_[xaxis, vdos[elem]]
            fname = elem + '_dos.pickle'
            with open(fname, 'wb') as wbfile:
                pickle.dump(data, wbfile)


if __name__ == "__main__":
    main()
