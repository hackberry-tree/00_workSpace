#!/opt/local/bin/python
# -*- coding: utf-8 -*-
"""
データfitting用のモジュール
"""
import numpy
import pylab
import scipy.optimize


class FitData(object):
    """
    fitting用のobject
    fitting_analysis.pyに移行
    """
    @classmethod
    def fit_arbfunc(cls, x, y, func, coefs):
        """
        任意の関数で fit
        memo として
        func の形式は func(coefs, x)
        1: 誤差の関数作成
        2: 最小二乗法
        """
        def err(coefs, x, y):
            return y - func(coefs, x)
        fit, fit_err = cls.leastsq(err, coefs, args=(x, y))
        return fit, fit_err


    @staticmethod
    def Stineman_interp_fit(x, y): #pylint: disable=C0103
        """
        データ補間曲線
        """
        yp = None #pylint: disable=C0103
        xi = pylab.linspace(x[0], x[-1], 100) #pylint: disable=E1101,C0103
        yi = pylab.stineman_interp(xi, x, y, yp) #pylint: disable=C0103
        return xi, yi

    @staticmethod
    def local_minimum(x, y): #pylint: disable=C0103
        """to find local minimum"""
        min_pos = (numpy.r_[True, y[1:] < y[:-1]] & #pylint: disable=E1101
                   numpy.r_[y[:-1] < y[1:], True]) #pylint: disable=E1101
        x_min = x[min_pos]
        y_min = y[min_pos]
        return x_min, y_min

    @classmethod
    def Murnaghan_fit(cls, volume_array, energy_array):
    #pylint: disable=C0111,C0103,R0914
        """
        Murnaghanの式でfitting
        vfit: 間隔を細かく取ったvolume
        cls.Murnaghan_func(): fitting結果のenergy
        murnumpyars
        """
        v = volume_array #pylint: disable=C0103
        e = energy_array #pylint: disable=C0103
        vfit = numpy.linspace(min(v), max(v), 100) #pylint: disable=E1101
        # fit a parabola to the data, y = ax^2 + bx + c
        a, b, c = pylab.polyfit(v, e, 2) #pylint: disable=E1101,C0103

        #now here are our initial guesses.
        v0 = -b/(2*a) #pylint: disable=C0103
        e0 = a*v0**2 + b*v0 + c #pylint: disable=C0103
        b0 = 2*a*v0 #pylint: disable=C0103
        bP = 4 #pylint: disable=C0103

        x0 = [e0, b0, bP, v0] #pylint: disable=C0103
        xfit, xerr = cls.leastsq(cls.Murnaghan_err, x0, args=(v, e))
        return vfit, cls.Murnaghan_func(xfit, vfit), xfit, xerr

    @staticmethod
    def Murnaghan_func(parameters, vol): #pylint: disable=C0103
        """
        given a vector of parameters and volumes, return a vector of energies.
        equation From PRB 28,5480 (1983)
        """
        E0 = parameters[0]
        B0 = parameters[1]
        BP = parameters[2]
        V0 = parameters[3]
        E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)
        return E

    @classmethod
    def Murnaghan_err(cls, pars, x, y):
        """
        Murnaghan fittingからのの誤差
        最小二乗法で使用
        """
        err = y - cls.Murnaghan_func(pars, x)
        return err

    @staticmethod
    def leastsq(func, x0, args):
        """
        scipyを使った最小二乗法
        誤差もreturnする
        誤差(func)を最小化する
        """
        xfit, xcov, infodict, errmsg, success = \
            scipy.optimize.leastsq(func, x0, args, full_output=1)
        s_sq = (func(xfit, *args) ** 2).sum() / len(args[0])
        try:
            xcov *= s_sq
            error = [numpy.absolute(xcov[i][i])**0.5 for i in range(len(x0))]
        except TypeError:
            error = [0 for i in range(len(x0))]
        return xfit, numpy.array(error)


