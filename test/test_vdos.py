#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for vdos"""
import os
import pylab
from scipy import hamming, hanning
from scipy.fftpack import fft, ifft
import unittest
from vdos import *
from liquid2phase import VDoS, get_s_con
__date__ = "Apr 3 2015"

TEST_PATH = ("/Users/enoki/Researches/Analysis/Codes/01_testRun/")

def test():
    """
    main
    """
    unittest.main()


class Test(unittest.TestCase):
    """
    テスト
    """
    path = os.path.join(TEST_PATH, 'vdos')

    def _test_make_poscars(self):
        path = os.path.join(self.path, 'Al100Cu100')
        atoms = Atoms.from_xdatcar(
            os.path.join(path, 'XDATCAR'), start=-300, end=-1)
        atoms.make_poscars(1, dst=path)


    def test_cn(self):
        path = os.path.join(self.path, 'Al140Fe60')
        atoms = Atoms.from_xdatcar(
            os.path.join(path, 'XDATCAR'), start=0, end=-1)
        # atoms.get_cn_single(0, 0, 0.22, -1, 0.001)
        n, frac = atoms.get_cn(0.25, 0.015)
        get_s_con(n, np.array([0.7, 0.3]), frac)



    def _test_vdos(self):
        path = os.path.join(self.path, 'Al100Cu100')
        atoms = Atoms.from_xdatcar(
            os.path.join(path, 'XDATCAR'), start=-300, end=-1)
        xaxis, vdos = atoms.get_vdos()
        for elem in vdos:
            pylab.plot(xaxis, vdos[elem])
        pylab.xlim(0, 30)
        pylab.show()

    def _test_vacf(self):
        """
        test
        """
        path = os.path.join(self.path, 'Cu200')
        data = []
        data2 = []
        atoms = Atoms.from_xdatcar(os.path.join(path, 'XDATCAR'),
                                   start=-300, end=-1)
        # atoms.save_pickle(os.path.join(self.path, 'vdos.pickle'))
        # pylab.plot(atoms.vdos['total'].real)
        pylab.plot(atoms.get_vacf().mean(axis=-1))
        pylab.show()

    def _test_vaverage(self):
        i = 1
        atoms = Atoms.from_xdatcar(os.path.join(self.path, 'XDATCAR'),
                                   start=i*500, end=(i+1)*500)
        pylab.plot(atoms.velocity[:,0])
        pylab.show()
        return
        print(np.min(atoms.velocity))
        print(np.max(atoms.velocity))
        for i in range(49):
            atoms = Atoms.from_xdatcar(os.path.join(self.path, 'XDATCAR'), start=i*1000, end=(i+1)*1000)
            print(atoms.get_averaged_speed())


    def _test_fft(self):
        """
        横軸の補正の仕方について考察
        横軸の補正は変換前のデータの横幅 L を使って 1/L 倍することで求める
        """
        data = np.arange(0, 1, 0.001)
        data_i = data / 0.001
        #print(data_i)
        ddata = 0.001
        print(len(np.fft.fft(np.cos(data * 4 * np.pi)).real))
        pylab.plot(np.fft.fft(np.cos(data * 4 * np.pi)))
        pylab.show()
        #print(2 * np.pi == data_i * 2 * np.pi * ddata )

    def _test_l2p(self):
        """
        liquid2phase.py と vdos.py を組み合わせて entropy とグラフを作成
        """
        path = os.path.join(self.path, 'Al140Fe60')
        atoms = Atoms.from_xdatcar(
            os.path.join(path, 'XDATCAR'),  start=-200, end=-1)
        volume = atoms.get_volume()
        pomas = {"Cu": 63.546, "Al": 26.982, "Fe": 55.847}
        print(volume)
        xaxis, vdos = atoms.get_vdos()
        ave_math = sum([pomas[x]*y
                        for x, y in zip(atoms.elements,
                                        atoms.num_atoms)]) / sum(atoms.num_atoms)
        total, total_geg, total_ges = 0, 0, 0
        for elem, num in zip(atoms.elements, atoms.num_atoms):
            print(elem)
            print(num)
            data = np.c_[xaxis, vdos[elem]]
            analysis = VDoS(data, 2300, num, volume, pomas[elem])
            # analysis = VDoS(data, 2300, 200, 3324, 63.546)

            gas = analysis.get_vdos_gas()(analysis.vdos[:, 0])
            solid = analysis.get_vdos_solid()
            pylab.plot(xaxis, vdos[elem])
            pylab.plot(xaxis, gas)
            pylab.plot(xaxis, solid)
            for i in range(len(xaxis)):
                print("{0} {1} {2} {3}".format(xaxis[i], vdos[elem][i], gas[i], solid[i]))
            correction = pomas[elem] * num / sum(atoms.num_atoms) / ave_math
            geg = analysis.get_entropy_gas() * correction
            ges = analysis.get_entropy_solid() * correction
            print(geg)
            print(ges)
            total_geg += geg
            total_ges += ges
            total += geg + ges
        print(total_geg)
        print(total_ges)
        print(total)
        pylab.xlim(0, 30)
        # pylab.show()
        pylab.savefig("/Users/enoki/output.eps")


if __name__ == '__main__':
    test()
