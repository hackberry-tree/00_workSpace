#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for parse_atat"""


import os
import unittest
from benchmarker import Benchmarker
from montecarlo import MonteCarlo, NaClXtal, FCCXtal
import numpy as np

__date__ = "Aug 4 2015"

TEST_PATH = ("/Users/enoki/Researches/Analysis/Codes/01_testRun/")


def main():
    """
    Test
    """
    unittest.main()


class Test(unittest.TestCase):  #pylint: disable=R0903
    """
    Test
    """
    PATH = os.path.join(TEST_PATH, 'montecarlo/')

    def test_get_entropy_TO(self):
        """
        de を記録して entropy を算出する
        """
        path = os.path.join(self.PATH, "AlCu/wien/TO")
        fcc = FCCXtal.from_pickle_ecis(os.path.join(path, 'cluster.pickle'),
                                       arrange='random', conc=0.25, size=30)
        mc = MonteCarlo(fcc, T=10000)
        mc.loop_fcc_micro_single(2000)


    def _test_flip_de(self):
        """
        flip de を total energy の変化量と比較
        一致すれば O.K.
        cell size が 2 以上ならうまくいく
        それより小さい場合は direct にエネルギーを求めた方が良い
        """
        path = os.path.join(self.PATH, "AlCu/voldep/4.0/")
        fcc = FCCXtal.from_pickle_ecis(os.path.join(path, 'cluster.pickle'),
                                       arrange='random', conc=0.8, size=2)

        # s=0, s=1 のサイトを無作為に抽出
        s0 = fcc.cell == 0
        select_s0 = np.random.choice(range((s0).sum()), 1, replace=False)
        site_s0 = (np.array(np.where(s0)))[:, select_s0[0]]

        s1 = fcc.cell == 1
        select_s1 = np.random.choice(
            range((s1).sum()), 1, replace=False)
        site_s1 = (np.array(np.where(s1)))[:, select_s1[0]]

        print("before energy")
        before = fcc.get_energy()
        print(before)
        pred_de = fcc.get_exchange_de_small(site_s0, site_s1)/fcc.size**3/4
        fcc.cell[tuple(site_s0.T)] = 1
        fcc.cell[tuple(site_s1.T)] = 0
        after = fcc.get_energy()
        print("after energy")
        print(after)
        print("de")
        de = after - before
        print(de)
        print("Predict de")
        print(pred_de)
        assert (de - pred_de) ** 2 < 1e-6

    def _test_mc_fcc_micoro_single(self):
        """
        loop_fcc_mciro_single の テスト
        benchmarker で速度計測
        """
        # path = os.path.join(self.PATH, "AlCu/voldep/4.0/")
        path = os.path.join(self.PATH, "AlCu/wien/")
        fcc = FCCXtal.from_pickle_ecis(os.path.join(path, 'cluster.pickle'),
                                       arrange='random', conc=0.90, size=10)
        print(len(fcc.ecis))
        mc = MonteCarlo(fcc, T=300)
        with Benchmarker(width=20, loop=3) as bench:
            @bench("simple")
            def _(bm): #pylint: disable=W0613,C0111
                mc.loop_fcc_micro_single(10)

    def _test_gs(self):
        """
        基底状態探索
        エネルギーのみを出力
        """
        path = os.path.join(self.PATH, "AlCu/voldep/4.0/")

        fcc = FCCXtal.from_pickle_ecis(os.path.join(path, 'cluster.pickle'),
                                       arrange='random', conc=0.90, size=10)
        fcc.from_pickle_cell(os.path.join(path, 'AlCu.pickle'))

        # mc = MonteCarlo(fcc, T=500)
        # mc.loop_fcc_micro_single(500)
        # 90 300
        mc = MonteCarlo(fcc, T=20)
        mc.loop_fcc_micro_single(10000)

        fcc.make_poscar(os.path.join(path, 'POSCAR'))
        fcc.save_cell(os.path.join(path, 'AlCu'))

    def _test_AlCu_mc_grand(self):
        """
        グランドカノニカルの計算テスト
        """
        fcc = FCCXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_AlCu_09.pickle"),
            arrange='random', conc=0.05, size=10)
        # fcc.from_pickle_cell(os.path.join(self.PATH, 'AlCu.pickle'))
        # eci_point = - 454.351346
        # -3250 @ 1500
        # -3800 @ 8000 0.92
        # -4000 @ 10000 0.92
        # -4700 @ 13500
        fcc.eci_point = 89.45
        mc = MonteCarlo(fcc, T=0)
        mc.loop_fcc_grand(20)

        fcc.make_poscar(os.path.join(self.PATH, 'POSCAR'))
        fcc.save_cell(os.path.join(self.PATH, 'AlCu'))

    def _test_AlCu_mc_grand_TO(self):
        """
        グランドカノニカルの計算テスト TO 近似
        """
        fcc = FCCXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_AlCu_TO.pickle"),
            arrange='random', conc=0.6, size=30)
        # fcc.from_pickle_cell(os.path.join(self.PATH, 'AlCu.pickle'))
        mc = MonteCarlo(fcc, T=1000)
        mc.loop_fcc_grand(100)

        fcc.make_poscar(os.path.join(self.PATH, 'POSCAR'))
        fcc.save_cell(os.path.join(self.PATH, 'AlCu'))

    def _test_flip_de_AlCutetra(self):
        """
        flip de を total energy の変化量と比較
        single サイト
        """
        fcc = FCCXtal.from_pickle_ecis(
            os.path.join(self.PATH, "AlCu/tetra_stress30/cluster.pickle"),
            arrange='FCC_A', size=6)
        fcc.cell[1, 2, 1, 3] = 0
        fcc.cell[1, 1, 1, 3] = 0
        fcc.cell[2, 2, 1, 0] = 0
        fcc.cell[1, 2, 1, 0] = 0

        pred_de = fcc.get_flip_de()[2, 2, 1, 3]/fcc.size**3/4

        print("before")
        before = fcc.get_energy()
        print(before)

        fcc.cell[2, 2, 1, 3] = 0
        # fcc.cell[1, 2, 1, 2] = 0


        print((fcc.get_flip_de() == np.min(fcc.get_flip_de())).sum())
        after = fcc.get_energy()
        print(after)
        print("de")
        de = after - before
        print(de)
        print("Predict de")
        print(pred_de)
        print((de - pred_de) ** 2 < 1e-6)

    def _test_order_AlCu_TO(self):
        """
        規則相のエネルギーを icvm と比較
        null clucter は含まれていないので注意
        """
        path = os.path.join(self.PATH, "AlCu/wien/TO")
        fcc = FCCXtal.from_pickle_ecis(os.path.join(path, 'cluster.pickle'),
                                       arrange='L12_A', size=10)
        print("L12_A3B: -40.7461")
        print(fcc.get_energy())

        fcc = FCCXtal.from_pickle_ecis(os.path.join(path, 'cluster.pickle'),
                                       arrange='L10', size=10)
        print("L10_AB -148.668")
        print(fcc.get_energy())

        fcc = FCCXtal.from_pickle_ecis(os.path.join(path, 'cluster.pickle'),
                                       arrange='L12_B', size=10)
        print("L12_AB3 -183.749")
        print(fcc.get_energy())

    def _test4(self):
        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int_Ti.pickle"),
            arrange='random', conc_int=0.001, conc_sub=0.999, size=15)
        # print((-nacl.get_flip_de_int()/80))

        mc = MonteCarlo(nacl, T=1000)

        # int -569.921196 sub -240.030646
        # nacl.sub_ecis[0] = -600.
        # nacl.int_ecis[0] = -500
        #mc.loop_flip_nacl(20)

        #mc.loop_pairflip_nacl(100)

        mc.loop_flip_nacl_fix_conc(100)

        nacl.make_poscar(os.path.join(self.PATH, 'POSCAR'))

    def _test_order_Cr(self):
        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int.pickle"),
            arrange='FCC_u', size=2)
        print('FCC Fe 57.3544')
        print(nacl.get_energy() + 171.363194)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int.pickle"),
            arrange='NaCl_uu', size=2)
        print('NaCl FeC 531.492')
        print(nacl.get_energy() + 171.363194)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int.pickle"),
            arrange='NaCl_du', size=2)
        print('NaCl CrC 115.585')
        print(nacl.get_energy() + 171.363194)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int.pickle"),
            arrange='FCC_d', size=2)
        print('FCC Cr 171.363')
        print(nacl.get_energy() + 171.363194)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int.pickle"),
            arrange='FCC_u', size=10)
        nacl.cell[0,0,0,0] = 0
        single_Cr = nacl.get_energy()
        nacl.cell[0,0,0,5] = 1
        pair_CrC = nacl.get_energy()
        nacl.cell[0,0,0,5] = 0
        nacl.cell[5,5,5,5] = 1
        solo_CrC = nacl.get_energy()
        de = pair_CrC - solo_CrC
        print(de * 10**3 * 4)
        print()

    def _test_order_Ti(self):
        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int_Ti.pickle"),
            arrange='FCC_u', size=2)
        print('FCC Fe 63.0657')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int_Ti.pickle"),
            arrange='NaCl_uu', size=2)
        print('NaCl FeC 532.114')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int_Ti.pickle"),
            arrange='NaCl_du', size=2)
        print('NaCl TiC -839.445')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int_Ti.pickle"),
            arrange='FCC_d', size=2)
        print('FCC Ti 54.8679')
        print(nacl.get_energy() + 54.8679061)
        print()

    def _test_clus_ene(self):
        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int_Ti.pickle"),
            arrange='FCC_u', size=5)
        nacl.cell[0,0,0,0] = 0
        single_Ti = nacl.get_energy()
        nacl.cell[0,0,0,5] = 1
        pair_TiC = nacl.get_energy()
        nacl.cell[0,0,0,5] = 0
        nacl.cell[3,3,3,5] = 1
        solo_TiC = nacl.get_energy()
        de = pair_TiC - solo_TiC
        print(de * 5**3 * 4)

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int_Ti_onlyNN.pickle"),
            arrange='FCC_u', size=5)
        nacl.cell[1,0,0,0] = 0
        single_Ti = nacl.get_energy()
        nacl.cell[1,0,0,5] = 1
        pair_TiC = nacl.get_energy()
        nacl.make_poscar(os.path.join(self.PATH, 'POSCAR'))
        nacl.cell[1,0,0,5] = 0
        nacl.cell[3,3,3,5] = 1
        solo_TiC = nacl.get_energy()
        de = pair_TiC - solo_TiC
        print(de * 5**3 * 4)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "cluster_int_Ti.pickle"),
            arrange='FCC_u', size=5)
        all_Fe = nacl.get_energy()
        nacl.cell[0,0,0,0] = 0
        single_Ti = nacl.get_energy()
        nacl.cell[0,0,0,5] = 1
        pair_TiC = nacl.get_energy()
        nacl.cell[0,0,0,0] = 1
        single_C = nacl.get_energy()
        de = pair_TiC - single_Ti - single_C + all_Fe
        print(de * 5**3 * 8 / 2)

if __name__ == '__main__':
    main()
