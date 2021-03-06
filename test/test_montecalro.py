#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for parse_atat"""


import os
import unittest
from benchmarker import Benchmarker
from montecarlo import MonteCarlo, NaClXtal, FCCXtal, CEMParser, SubLattXtal
from montecarlo import QuadSite, NaClSite, BcciOctaSite, CEMParser02
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

    def _test_get_cluser_energy(self):
        """
        全てのクラスターを作成分解してエネルギーを出す
        """
        path = os.path.join(self.PATH, "i-s/fcci/Ti1C1/R6N4")
        # path = os.path.join(self.PATH, "i-s/fcci/Ti0C1/R4N4")
        # path = os.path.join(self.PATH, "i-s/fcci/Ti1C1/R4N4")
        clus = CEMParser.parse_logtxt_2sub(os.path.join(path, 'log.txt'))
        print(clus[0])
        clus_idxs = (
            np.array([[NaClSite.conv_site2idx(site) for site in cluster]
                for cluster in clus[0]]))
        fcci_xtal = NaClXtal.from_pickle_ecis(os.path.join(path, 'cluster.pickle'),
                                             size=5)
        # 全ての POSCAR を作成
        # for i in range(len(clus_idxs)):
        #     fcci_xtal.cell[:, :, :, :] = 0
        #     fcci_xtal.cell[tuple(np.array(clus_idxs[i]).T)] = 1
        #     fcci_xtal.make_poscar(os.path.join(path, 'POSCAR_{0}'.format(i)), 1, 1)

        # 全てのエネルギーを計算
        # fcci_xtal.sub_ecis[3] = 0
        # fcci_xtal.int_ecis[3] = 0
        # fcci_xtal.sub_ecis[4] = 0
        print(fcci_xtal.int_ecis)

        c_spin = 1
        a_spin = 1
        fcci_xtal.cell[:, :, :, :] = 0
        fcci_xtal.cell[:, :, :, 0:4] = 1 - c_spin
        ref_null = fcci_xtal.get_energy() * fcci_xtal.size ** 3 * 8
        print(ref_null )

        fcci_xtal.cell[:, :, :, :] = 0
        fcci_xtal.cell[tuple(np.array(clus_idxs[0]).T)] = 1
        fcci_xtal.cell[:, :, :, 0:4] = \
            (1 - fcci_xtal.cell[:, :, :, 0:4]) * (1 - c_spin) + \
            fcci_xtal.cell[:, :, :, 0:4] * c_spin
        ref_i = fcci_xtal.get_energy() * fcci_xtal.size ** 3 * 8 - ref_null

        fcci_xtal.cell[:, :, :, :] = 0
        fcci_xtal.cell[tuple(np.array(clus_idxs[1]).T)] = 1
        fcci_xtal.cell[:, :, :, 0:4] = \
            (1 - fcci_xtal.cell[:, :, :, 0:4]) * (1 - c_spin) + \
            fcci_xtal.cell[:, :, :, 0:4] * c_spin
        ref_s = fcci_xtal.get_energy() * fcci_xtal.size ** 3 * 8 - ref_null

        for i in range(len(clus_idxs)):
            fcci_xtal.cell[:, :, :, :] = 0
            fcci_xtal.cell[tuple(np.array(clus_idxs[i]).T)] = 1
            fcci_xtal.cell[:, :, :, 0:4] = \
                (1 - fcci_xtal.cell[:, :, :, 0:4]) * (1 - c_spin) + \
                fcci_xtal.cell[:, :, :, 0:4] * c_spin
            fcci_xtal.make_poscar(os.path.join(path, 'POSCAR_{0}'.format(i+1)), 1, 1)
            fcci_xtal.cell[:, :, :, :] = c_spin
            fcci_xtal.cell[:, :, :, 4:] = 1
            tot_e = fcci_xtal.get_energy() * fcci_xtal.size ** 3 * 8
            num_s = ((1 - fcci_xtal.get_conc()['subs']) * (1 - c_spin) +
                     fcci_xtal.get_conc()['subs'] * c_spin) * fcci_xtal.size ** 3 * 4
            num_i = fcci_xtal.get_conc()['inter'] * fcci_xtal.size ** 3 * 4
            num = (num_s + num_i)
            clus_e = (tot_e - num_s * ref_s - num_i * ref_i - ref_null) / num
            print(i+1, num_s, num_i, clus_e)
            # clus_e2 = (tot_e - )
        fcci_xtal.cell[:, :, :, :] = 1

    def _test_mc_bcci2(self):
        path = os.path.join(self.PATH, "i-s/bcci/TiC_R6N4_330str")
        CEMParser.from_dirc_bcci(path, BcciOctaSite)

    def _test_mc_bcci(self):
        path = os.path.join(self.PATH, "i-s/bcci/TiC_spin_r_R6N4")
        path = os.path.join(self.PATH, "i-s/bcci/TiC_R6N4_330str")
        CEMParser02.from_dirc_bcci(path, BcciOctaSite)
        CrC = SubLattXtal.from_pickle_ecis(
            5, os.path.join(path, 'cluster.pickle'),
            BcciOctaSite, arrange='random', conc=(0.1, 0.1))
        # CrC = SubLattXtal(2, out, BcciOctaSite, arrange='A3C', conc=(0.1, 0.1))

        CrC.make_poscar(os.path.join(self.PATH, "i-s/bcci/TiC_R6N4_330str/POSCAR_init"), 1, 1)

        mc = MonteCarlo(CrC, T=1000)
        e = mc.loop_bcci_micro_single(1)
        CrC.make_poscar(os.path.join(self.PATH, "i-s/bcci/TiC_R6N4_330str/POSCAR"), 1, 1)
        print(e)

    def _test_load_poscar(self):
        """POSCAR から cell を作成"""
        path = os.path.join(self.PATH, "i-s/poscar")
        # fcci = NaClXtal.from_pickle_ecis(os.path.join(path, "cluster.pickle"),
        #                                   arrange='random', size=10, conc=(0.01, 0.01))
        # fcci.make_poscar(os.path.join(path, "POSCAR"), 1, 1)
        fcci = NaClXtal.from_pickle_ecis(os.path.join(path, "cluster.pickle"),
                                         arrange='FCC_d', size=10)

        fcci.load_poscar(os.path.join(path, "POSCAR"), 1, 1)
        fcci.make_poscar(os.path.join(path, "POSCAR_rep"), 1, 1)

        fcci = NaClXtal.from_pickle_ecis(os.path.join(path, "cluster.pickle"),
                                         arrange='FCC_d', size=20)
        fcci.load_poscar(os.path.join(path, "POSCAR_4840"), 1, 1)
        fcci.make_poscar(os.path.join(path, "POSCAR_rep2"), 0, 1, True)

    def _test_compare_pos_vs_clus_bcci(self):
        """
        bcciにて
        log.txtから与えられた元素種と site 位置から判別した元素種とが一致するか
        check 用の method
        """
        path = os.path.join(self.PATH, "i-s/bcci/CrC")
        judges = CEMParser.compare_pos_vs_clus(
            os.path.join(path, 'log.txt'), BcciOctaSite)
        for judge in judges:
            assert judge

    def _test_bcci(self):
        path = os.path.join(self.PATH, "i-s/bcci/CrC")
        # clus = CEMParser.parse_logtxt_bcci(os.path.join(path, 'log.txt'))
        # # print(clus[2])
        # # CEMParser.symm_cubic_bcci(clus[0], BcciOctaSite)
        # clusters = []
        # for c in clus:
        #     clusters.append(CEMParser.symm_cubic_bcci(c, BcciOctaSite))
        # ecis = CEMParser.parse_ecitxt(os.path.join(path, 'eci.txt'))
        # out = {}
        # for label in ['C', 'A1', 'A2', 'A3']:
        #     tmp_clus = []
        #     tmp_ecis = []
        #     for i in sorted(ecis.keys()):
        #         if clusters[i][label]:
        #             tmp_clus.append(
        #                 np.array([[BcciOctaSite.conv_site2idx(site)
        #                            for site in cluster]
        #                           for cluster in clusters[i][label]]))
        #             tmp_ecis.append(ecis[i])
        #     out.update({label:[tmp_clus, tmp_ecis]})
        # print(out['C'])
        # CEMParser02.from_dirc_bcci(path, BcciOctaSite)
        # CrC = SubLattXtal.from_pickle_ecis(
        #     5, os.path.join(path, 'cluster.pickle'),
        #     BcciOctaSite, arrange='random', conc=(0.1, 0.1))
        pickle = os.path.join(path, 'cluster.pickle')

        print("A3C -.678871E+04")
        print(SubLattXtal.from_pickle_ecis(
            1, pickle, BcciOctaSite, arrange='A3C').get_energy() - 9614.98293)

        print("A2C -.685326E+04")
        print(SubLattXtal.from_pickle_ecis(
            1, pickle, BcciOctaSite, arrange='A2C').get_energy() - 9614.98293)

        print("A1C -.734503E+04")
        print(SubLattXtal.from_pickle_ecis(
            1, pickle, BcciOctaSite, arrange='A1C').get_energy() - 9614.98293)

        print("BCC_C -.829233E+04")
        print(SubLattXtal.from_pickle_ecis(
            1, pickle, BcciOctaSite, arrange='BCC_u').get_energy() - 9614.98293)

        print("A3D -.705443E+04")
        print(SubLattXtal.from_pickle_ecis(
            1, pickle, BcciOctaSite, arrange='A3D').get_energy() - 9614.98293)

        print("A2D -.723986E+04")
        print(SubLattXtal.from_pickle_ecis(
            1, pickle, BcciOctaSite, arrange='A2D').get_energy() - 9614.98293)

        print("A1D -.832914E+04")
        print(SubLattXtal.from_pickle_ecis(
            1, pickle, BcciOctaSite, arrange='A1D').get_energy() - 9614.98293)

        print("BCC_D -.961498E+04")
        print(SubLattXtal.from_pickle_ecis(
            1, pickle, BcciOctaSite, arrange='BCC_d').get_energy() - 9614.98293)

        # print("ACD")
        # print(SubLattXtal.from_pickle_ecis(
        #    4, pickle, BcciOctaSite, arrange='ACD').get_energy() - 9614.98293)


        CrC = SubLattXtal.from_pickle_ecis(
            4, pickle, BcciOctaSite, arrange='BCC_u')
        CrC.cell[0,0,0,2] = 1
        CrC.cell[0,0,0,0] = 0
        crc1 = CrC.get_energy()
        CrC.cell[0,0,0,2] = 0
        CrC.cell[2,2,2,2] = 1
        crc2 = CrC.get_energy()
        print((crc1 - crc2)*4**3*8)

    def test_exch_de_bcci(self):
        path = os.path.join(self.PATH, "i-s/bcci/CrC")
        pickle = os.path.join(path, 'cluster.pickle')

        print("A3C -.678871E+04")
        print(SubLattXtal.from_pickle_ecis(
            1, pickle, BcciOctaSite, arrange='A3C').get_energy() - 9614.98293)

        bcci = SubLattXtal.from_pickle_ecis(
             30, pickle, BcciOctaSite, arrange='random', conc=(0.1,0.1))
        bcci.cell[0,0,0,7] = 1
        # print(bcci.get_energy())
        print(bcci.get_energy())
        with Benchmarker(width=20, loop=10) as bench:
            @bench("get_energy")
            def _(bm): #pylint: disable=W0613,C0111
                bcci.get_exchange_de(np.array([0,0,0,0]), np.array([0,0,0,0]))
            @bench("get_energy2")
            def __(bm): #pylint: disable=W0613,C0111
                bcci.get_exchange_de2(np.array([0,0,0,0]), np.array([0,0,0,0]))

        return
        for i in range(10):
            s_or_i = np.random.choice([0, 1], 1)[0]
            print(['sub', 'int'][s_or_i])

            s0 = [bcci.cell[:, :, :, 0:2] == 0,
                  bcci.cell[:, :, :, 2:] == 0][s_or_i]
            select_s0 = np.random.choice(range((s0).sum()), 1, replace=False)
            site_s0 = ((np.array(np.where(s0)))[:, select_s0[0]] +
                       np.array([0,0,0,s_or_i*2]))

            s1 = [bcci.cell[:, :, :, 0:2] == 1,
                  bcci.cell[:, :, :, 2:] == 1][s_or_i]
            select_s1 = np.random.choice(range((s1).sum()), 1, replace=False)
            site_s1 = ((np.array(np.where(s1)))[:, select_s1[0]] +
                       np.array([0,0,0,s_or_i*2]))
            de = bcci.get_exchange_de2(site_s0, site_s1)
            e0 = bcci.get_energy()
            bcci.cell[tuple(site_s0.T)] = 1
            bcci.cell[tuple(site_s1.T)] = 0
            e1 = bcci.get_energy()
            print((de - (e1-e0)*bcci.size**3*8) ** 2 < 1e-16)

        # bcci = SubLattXtal.from_pickle_ecis           (
        #     3, pickle, BcciOctaSite, arrange='random')
        # while bcci.cell[0,0,0,0] != 0 or bcci.cell[0,0,0,1] != 1:
        #     bcci = SubLattXtal.from_pickle_ecis(
        #         3, pickle, BcciOctaSite, arrange='random')
        # print(bcci.cell[0,0,0,0])
        # print(bcci.cell[0,0,0,1])
        # e0 = bcci.get_energy()
        # de = bcci.get_exchange_de2(np.array([0,0,0,0]),
        #                            np.array([0,0,0,1]))
        # bcci.cell[0,0,0,0] = 1
        # bcci.cell[0,0,0,1] = 0
        # e1 = bcci.get_energy()
        # print(de)
        # print((e1-e0)*bcci.size**3*8)

        # while bcci.cell[0,0,0,2] != 0 or bcci.cell[0,0,0,5] != 1:
        #     bcci = SubLattXtal.from_pickle_ecis(
        #         3, pickle, BcciOctaSite, arrange='random')
        # print(bcci.cell[0,0,0,2])
        # print(bcci.cell[0,0,0,5])
        # e0 = bcci.get_energy()
        # de = bcci.get_exchange_de2(np.array([0,0,0,2]),
        #                            np.array([0,0,0,5]))
        # bcci.cell[0,0,0,2] = 1
        # bcci.cell[0,0,0,5] = 0
        # e1 = bcci.get_energy()
        # print(de)
        # print((e1-e0)*bcci.size**3*8)

        # while bcci.cell[0,0,0,2] != 0 or bcci.cell[0,0,0,3] != 1:
        #     bcci = SubLattXtal.from_pickle_ecis(
        #         3, pickle, BcciOctaSite, arrange='random')
        # print(bcci.cell[0,0,0,2])
        # print(bcci.cell[0,0,0,3])
        # e0 = bcci.get_energy()
        # de = bcci.get_exchange_de2(np.array([0,0,0,2]),
        #                            np.array([0,0,0,3]))
        # bcci.cell[0,0,0,2] = 1
        # bcci.cell[0,0,0,3] = 0
        # e1 = bcci.get_energy()
        # print(de)
        # print((e1-e0)*bcci.size**3*8)

    def _test_bcci2(self):
        path = os.path.join(self.PATH, "i-s/bcci/CrC_spin_r")
        clus = CEMParser.parse_logtxt_bcci(os.path.join(path, 'log.txt'))
        # print(clus[2])
        # CEMParser.symm_cubic_bcci(clus[0], BcciOctaSite)
        clusters = []
        for c in clus:
            clusters.append(CEMParser.symm_cubic_bcci(c, BcciOctaSite))
        ecis = CEMParser.parse_ecitxt(os.path.join(path, 'eci.txt'))
        out = {}
        for label in ['C', 'A1', 'A2', 'A3']:
            tmp_clus = []
            tmp_ecis = []
            for i in sorted(ecis.keys()):
                if clusters[i][label]:
                    tmp_clus.append(
                        np.array([[BcciOctaSite.conv_site2idex(site)
                                   for site in cluster]
                                  for cluster in clusters[i][label]]))
                    tmp_ecis.append(ecis[i])
            out.update({label:[tmp_clus, tmp_ecis]})
        # print(out['C'])

        print(SubLattXtal(1, out, BcciOctaSite, arrange='A3C').get_energy() - 8292.331923)
        print(SubLattXtal(1, out, BcciOctaSite, arrange='A2C').get_energy() - 8292.33192)
        print(SubLattXtal(1, out, BcciOctaSite, arrange='A1C').get_energy() - 8292.33192)
        print(SubLattXtal(1, out, BcciOctaSite, arrange='BCC_u').get_energy() - 8292.33192)
        print(SubLattXtal(1, out, BcciOctaSite, arrange='A3D').get_energy() - 8292.33192)
        print(SubLattXtal(1, out, BcciOctaSite, arrange='A2D').get_energy() - 8292.33192)

        print(SubLattXtal(1, out, BcciOctaSite, arrange='A1D').get_energy() - 8292.33192)

        print(SubLattXtal(1, out, BcciOctaSite, arrange='BCC_d').get_energy() - 8292.33192)

        print(SubLattXtal(4, out, BcciOctaSite, arrange='ACD').get_energy() - 8292.33192)

        CrC = SubLattXtal(4, out, BcciOctaSite, arrange='BCC_d')
        CrC.cell[0,0,0,2] = 1
        CrC.cell[0,0,0,0] = 1
        crc1 = CrC.get_energy()
        CrC.cell[0,0,0,2] = 0
        CrC.cell[2,2,2,2] = 1
        crc2 = CrC.get_energy()
        print((crc1 - crc2)*4**3*8)

    def _test_bcci3(self):
        path = os.path.join(self.PATH, "i-s/bcci/TiC_spin_r_R6N4")
        path = os.path.join(self.PATH, "i-s/bcci/TiC_R6N4_330str")
        clus = CEMParser.parse_logtxt_bcci(os.path.join(path, 'log.txt'))
        # print(clus[2])
        # CEMParser.symm_cubic_bcci(clus[0], BcciOctaSite)
        clusters = []
        for c in clus:
            clusters.append(CEMParser.symm_cubic_bcci(c, BcciOctaSite))
        ecis = CEMParser.parse_ecitxt(os.path.join(path, 'eci.txt'))
        out = {}
        for label in ['C', 'A1', 'A2', 'A3']:
            tmp_clus = []
            tmp_ecis = []
            for i in sorted(ecis.keys()):
                if clusters[i][label]:
                    tmp_clus.append(
                        np.array([[BcciOctaSite.conv_site2idex(site)
                                   for site in cluster]
                                  for cluster in clusters[i][label]]))
                    tmp_ecis.append(ecis[i])
            out.update({label:[tmp_clus, tmp_ecis]})
        # print(out['C'])

        print(SubLattXtal(1, out, BcciOctaSite, arrange='A3C').get_energy() - 7962.21404)
        print(SubLattXtal(1, out, BcciOctaSite, arrange='A2C').get_energy() - 7962.21404)
        print(SubLattXtal(1, out, BcciOctaSite, arrange='A1C').get_energy() - 7962.21404)
        print(SubLattXtal(1, out, BcciOctaSite, arrange='BCC_u').get_energy() - 7962.21404)
        print(SubLattXtal(1, out, BcciOctaSite, arrange='A3D').get_energy() - 7962.21404)
        print(SubLattXtal(1, out, BcciOctaSite, arrange='A2D').get_energy() - 7962.21404)

        print(SubLattXtal(1, out, BcciOctaSite, arrange='A1D').get_energy() - 7962.21404)

        print(SubLattXtal(1, out, BcciOctaSite, arrange='BCC_d').get_energy() - 7962.21404)

        print(SubLattXtal(4, out, BcciOctaSite, arrange='ACD').get_energy() - 7962.21404)

        CrC = SubLattXtal(4, out, BcciOctaSite, arrange='BCC_u')
        CrC.cell[0,0,0,2] = 1
        CrC.cell[0,0,0,0] = 0
        crc1 = CrC.get_energy()
        CrC.cell[0,0,0,2] = 0
        CrC.cell[2,2,2,7] = 1
        crc2 = CrC.get_energy()
        print((crc1 - crc2)*4**3*8)

    def _test_order_Ti(self):
        """
        Fe-Ti-C の test
        """
        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti0C1/cluster.pickle"),
            arrange='FCC_u', size=2)
        print('FCC Fe: 63.0657')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti0C1/cluster.pickle"),
            arrange='NaCl_uu', size=2)
        print('NaCl FeC: 532.114')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti0C1/cluster.pickle"),
            arrange='NaCl_du', size=2)
        print('NaCl TiC: -839.445')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti0C1/cluster.pickle"),
            arrange='FCC_d', size=2)
        print('FCC Ti: 54.8679')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti0C1/cluster.pickle"),
            arrange='order_00', size=2)
        print('FCC Ti: -243.445')
        print(nacl.get_energy() + 54.8679061)
        nacl.make_poscar(os.path.join(self.PATH, "i-s/fcci/Ti0C1/POSCAR_test"), 0, 1, True)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti0C1/cluster.pickle"),
            arrange='order_01', size=2)
        print('FCC Ti: -175.298')
        print(nacl.get_energy() + 54.8679061)
        nacl.make_poscar(os.path.join(self.PATH, "i-s/fcci/Ti0C1/POSCAR_test"), 0, 1, True)
        print()

    def _test_order_Ti_rev(self):
        """
        Fe-Ti-C の test
        """
        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti1C1/cluster.pickle"),
            arrange='FCC_d', size=2)
        print('FCC Fe: 63.0657')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti1C1/cluster.pickle"),
            arrange='NaCl_du', size=2)
        print('NaCl FeC: 532.114')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti1C1/cluster.pickle"),
            arrange='NaCl_uu', size=2)
        print('NaCl TiC: -839.445')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti1C1/cluster.pickle"),
            arrange='FCC_u', size=2)
        print('FCC Ti: 54.8679')
        print(nacl.get_energy() + 54.8679061)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti1C1/cluster.pickle"),
            arrange='order_01', size=2)
        print('FCC Ti: -243.445')
        print(nacl.get_energy() + 54.8679061)
        nacl.make_poscar(os.path.join(self.PATH, "i-s/fcci/Ti1C1/POSCAR_test"), 0, 1, True)
        print()

        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti1C1/cluster.pickle"),
            arrange='order_00', size=2)
        print('FCC Ti: -175.298')
        print(nacl.get_energy() + 54.8679061)
        nacl.make_poscar(os.path.join(self.PATH, "i-s/fcci/Ti1C1/POSCAR_test"), 0, 1, True)
        print()

    def _test_mc_micro_nacl(self):
        """
        i-s 系での mc 計算のてすと
        """
        nacl = NaClXtal.from_pickle_ecis(
            os.path.join(self.PATH, "i-s/fcci/Ti1C1/cluster.pickle"),
            arrange='random', conc=(0.001, 0.999), size=2)
        # print((-nacl.get_flip_de_int()/80))

        mc = MonteCarlo(nacl, T=1000)
        # mc.loop_flip_nacl_fix_conc(100)
        e = mc.loop_nacl_micro(10)
        lines = ""
        for i in range(len(e[0])):
            lines += str(i) + "\t" + str(e[0][i][2]) + "\t"
            for val in e[1][i][0]:
                lines += str(val) + "\t"
            for val in e[1][i][1]:
                lines += str(val) + "\t"
            lines += "\n"
        print(lines)


        # nacl.make_poscar(os.path.join(self.PATH, 'POSCAR'))

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

    def _test_CEM_Parser(self):
        """
        CEM Parser をテストする
        """
        path = "/Users/enoki/Researches/analysis/codes/01_testRun/montecarlo/AlCu/warren/150910_recalc"
        CEMParser.from_dirc_std(path, QuadSite)

        path = os.path.join(self.PATH, "AlCu/wien/TO")
        CEMParser.from_dirc_std(path, QuadSite)
        path = os.path.join(self.PATH, "AlCu/tetra_2R2N")
        CEMParser.from_dirc_prim_fcc(path, QuadSite)
        path = os.path.join(self.PATH, "i-s/fcci/Ti0C1/R4N4")
        CEMParser.from_dirc_2sub(path, NaClSite)
        # path = os.path.join(self.PATH, "i-s/bcci/CrC")
        # CEMParser.from_dirc_2sub(path, NaClSite)
        path = os.path.join(self.PATH, "i-s/bcci/TiC_R6N4_330str")
        CEMParser.from_dirc_2sub(path, NaClSite)

    def _test_get_entropy_TO(self):
        """
        de を記録して entropy を算出する
        """
        path = os.path.join(self.PATH, "AlCu/wien/TO")
        fcc = FCCXtal.from_pickle_ecis(os.path.join(path, 'cluster.pickle'),
                                       arrange='random', conc=0.25, size=30)
        mc = MonteCarlo(fcc, T=10000)
        mc.loop_fcc_micro_single(2000)

    def _test_flip_de_fcc(self):
        """
        flip de を total energy の変化量と比較
        一致すれば O.K.
        cell size が 3 以上ならうまくいく
        それより小さい場合は direct にエネルギーを求めた方が良い
        """
        path = os.path.join(self.PATH, "AlCu/voldep/4.0/")
        # path = os.path.join(self.PATH, "AlCu/wien/TO/")
        fcc = FCCXtal.from_pickle_ecis(os.path.join(path, 'cluster.pickle'),
                                       arrange='random', conc=0.8, size=3)
        for _ in range(10):
            # s=0, s=1 のサイトを無作為に抽出
            print("random flip")
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
            pred_de = fcc.get_exchange_de(site_s0, site_s1)/fcc.size**3/4
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
