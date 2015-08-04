#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parameters from
K. Terakura et al. PRB, 35, 2169
"""
from __future__ import division
import pylab

class FittingParam(object):
    """
    fitting parametersを格納
    """
    def __init__(self, p, q, r, a0, dE, conc):
        self.p = p
        self.q = q
        self.r = r
        self.a0 = a0
        self.dE = dE
        self.conc = conc

        self.p2n = p ** 7
        self.qn = q ** 3.5

    def energy(self, a):
        """
        格子長aにおけるエネルギーをreturn
        """
        a1 = a ** 7
        a2 = a ** 3.5
        return self.p2n / a1 - self.qn / a2 + self.r

class ECI(object):
    """
    ECIを格納
    """
    def __init__(self, fitting_params):
        self.fitparams = {}
        for param in fitting_params:
            self.fitparams.update({param.conc: param})

    def _calc_v(self, a, norm, c0, c1, c2, c3, c4):
        v = c0 * self.fitparams[0/4.].energy(a)
        v += c1 * self.fitparams[1/4.].energy(a)
        v += c2 * self.fitparams[2/4.].energy(a)
        v += c3 * self.fitparams[3/4.].energy(a)
        v += c4 * self.fitparams[4/4.].energy(a)
        v /= norm
        return v

    def v0(self, a):
        return self._calc_v(a, 16, 1, 4, 6, 4, 1)

    def v1(self, a):
        return self._calc_v(a, 4, 1, 2, 0, -2, -1)

    def v2(self, a):
        return self._calc_v(a, 8, 3, 0, -6, 0, 3)

    def v3(self, a):
        return self._calc_v(a, 4, 1, -2, 0, 2, -1)

    def v4(self, a):
        return self._calc_v(a, 16, 1, -4, 6, -4, 1)

class P_Sawada(object):
    """
    sawadaさんのプログラムのパラメータ
    3*3のcの値を入力
    """
    def __init__(self, c):
        self.c = c

    def v(self, j, a):
        a1 = a ** 3.5
        a2 = a ** 7
        v = self.c[j][0]
        v += self.c[j][1] / a1
        v += self.c[j][2] / a2
        # v /= 8
        v *= 1000
        return v




def main():
    param = {}
    Au = FittingParam(8.584, 11.677, 2.155, 7.692, 0.0, 0)
    print(Au.energy(Au.a0))
    CuAu3 = FittingParam(8.217, 11.013, 1.932, 7.474, -0.0100, 1/4.)
    print(CuAu3.energy(CuAu3.a0))
    Cu2Au2 = FittingParam(7.877, 10.450, 1.788, 7.237, -0.0205, 0.5)
    print(Cu2Au2.energy(Cu2Au2.a0))
    Cu3Au = FittingParam(7.474, 9.749, 1.586, 6.986, -0.0191, 3/4.)
    print(Cu3Au.energy(Cu3Au.a0))
    Cu = FittingParam(7.055, 9.033, 1.410, 6.717, 0.0, 1)
    print(Cu.energy(Cu.a0))
    params = [Au, CuAu3, Cu2Au2, Cu3Au, Cu]
    eci = ECI(params)
    print(eci.v0(7.8))

    list_a = [x / 100 for x in range(660, 780)]
    list_v0 = [eci.v0(x) * 1000 for x in list_a]
    list_v1 = [eci.v1(x) * 1000 for x in list_a]
    list_v2 = [eci.v2(x) * 1000 for x in list_a]
    list_v3 = [eci.v3(x) * 1000 for x in list_a]
    list_v4 = [eci.v4(x) * 1000 for x in list_a]

    c = []
    c.append([1.77281, -3693.31368, 1932616.12707])
    c.append([0.35925, -1576.20925, 1254312.59500])
    c.append([-0.00413, -104.24365, 202856.45586])
    c.append([0.01325, -36.57425, 27836.99500])
    c.append([0.01381, -30.42018, 16596.62707])
    vs = P_Sawada(c)
    list_v0s = [vs.v(0, x) for x in list_a]
    list_v1s = [vs.v(1, x) for x in list_a]
    list_v2s = [vs.v(2, x) for x in list_a]
    list_v3s = [vs.v(3, x) for x in list_a]
    list_v4s = [vs.v(4, x) for x in list_a]

    pylab.plot(list_a, list_v2)
    pylab.plot(list_a, list_v2s, 'x')

    pylab.show()


if __name__ == '__main__':
    main()
