#!/usr/bin/env python
"""
積分計算をtestする
"""
import os
import math
from scipy.integrate import quad
import numpy as np
import pylab


PATH = "/Users/enoki/Desktop/"

fname = os.path.join(PATH, 'VDOS.txt')
with open(fname, 'r') as rfile:
    read_lines = rfile.readlines()

data1 = []
for line in read_lines[1:]:
    data1.append([float(x) for x in line.split()])

data = np.array(data1)

#print(data)

def senkei_hokan(a, b):
    """
    a = [xa, ya]
    b = [xb, yb]
    が与えられたときその間を線形で結んだ直線f(x)をreturnする
    f(xa) = ya, f(xb) = yb
    """
    _a = np.array(a)
    _b = np.array(b)
    delta = _b - _a
    alpha = delta[1] / delta[0]
    closs = _a * _b[::-1]
    beta =  (closs[1] - closs[0]) / delta[0]
    return lambda x: alpha * x + beta

total = 0
err = 0
for i in range(len(data)-1):
    y, abserr = quad(
        senkei_hokan(data[i][0:2], data[i+1][0:2]), data[i][0], data[i+1][0])
    total += y
    err += abserr

def delta(T, N, V, M, D):
    """
    融解の無次化された係数Δを算出
    T: 温度 (K)
    N: 原子数
    V: 体積 (Å^3)
    M: 質量 (原子量)
    D: v=0におけるVDOS (THz^-1)
    """
    NA = 6.02e23
    _M = M / NA / 1000 # 単位kg
    kb = 1.3806488e-23
    beta = 1 / kb / T # 単位J
    _D = D * 10 ** -12 # 単位 /s
    _V = V * 10 ** -30
    return 2 * _D / 9 / N * ((math.pi / _M / beta) ** (0.5)) * \
        ((N / _V) ** (1. / 3)) * ((6 / math.pi) ** (2. / 3))

#print(delta(2300, 200, 4311.57, 26.9815395, 89.8881))
def dos_gas(d0, F, N):
    def vdos(v):
        return d0 / (1 + ((math.pi * d0 * v) / ( 6 * F * N)) ** 2 )
    return vdos

def dos_sol(DOS, dos_g):
    return DOS[:, 1] - dos_g(DOS[:, 0])

def weight_solid(T):
    kb = 1.3806488e-23
    beta = 1 / kb / T
    h = 6.62606957e-34# 単位Js
    def ws(v):
        bhv = beta * h * v * 10 ** 12
        return bhv / (np.exp(bhv) - 1) - np.log(1 - np.exp(-bhv))
    return ws

def packing_fraction(F, delta):
    return (F / delta) ** (2/3)

y = packing_fraction(0.501, 0.821)

f = np.linspace(7.5e-3, 30, num=50)

ws = weight_solid(2300)

pylab.plot(f, ws(f))
#pylab.show()




gas = dos_gas(90, 0.501, 200)
dgas = gas(f)

dsol = dos_sol(data, gas)
times = dsol * ws(data[:, 0])
print(ws(data[:, 0]))
print(times)

# S_solを算出
total = 0
for i in range(len(data)-1):
    y, abserr = quad(
        senkei_hokan([data[i][0], times[i]], [data[i+1][0], times[i+1]]), data[i][0], data[i+1][0])
    total += y
    err += abserr
print(total*8.31/200)

# 理想気体のエントロピーs_ig
def s_ig(N, V, M, T):
    kb = 1.3806488e-23
    beta = 1 / kb / T
    NA = 6.02e23
    _M = M / NA / 1000 # 単位kg
    _V = V * 10 ** -30
    h = 6.62606957e-34# 単位Js
    cT = (2 * math.pi * _M / beta / h / h) ** (3. / 2)
    return (N * kb * np.log(_V * cT / N) + 5. / 2 * N * kb) / N

def s_ig2(N, V, M, T):
    kb = 1.3806488e-23
    beta = 1 / kb / T
    NA = 6.02e23
    _M = M / NA / 1000 # 単位kg
    _V = V * 10 ** -30
    h = 6.62606957e-34# 単位Js
    def s(v):
        e = h * v * 1e15
        return kb * N * (3. / 2 * np.log(4 * math.pi * _M * e / 3 / N / h / h) \
            + np.log(_V/N) + 5 / 2) / N
    return s



def weight_gas(sig, F, delta):
    y = (F / delta) ** (2/3)
    fy = F * y
    kb = 1.3806488e-23
    second = np.log((1 + fy + (fy) ** 2 - (fy) ** 3) / (1 - fy) ** 3)
    third = (fy * (3 * fy - 4)) / ((1 - fy) ** 2)
    return 1./3 * (sig/kb + second + third)

def weight_gas2(sig, F, delta):
    y = (F / delta) ** (2/3)
    fy = F * y
    kb = 1.3806488e-23
    second = np.log((1 + fy + (fy) ** 2 - (fy) ** 3) / (1 - fy) ** 3)
    third = (fy * (3 * fy - 4)) / ((1 - fy) ** 2)
    def s(v):
        return 1./3 * (sig(v)/kb + second + third)
    return s

sig = s_ig(200, 4311.57, 26.9815395, 2300)
print(sig/1.3806488e-23)
wg = weight_gas(sig, 0.501, 0.821)
times2 = wg * gas(data[:, 0])
total = 0
for i in range(len(data)-1):
    y, abserr = quad(
        senkei_hokan([data[i][0], times2[i]], [data[i+1][0], times2[i+1]]), data[i][0], data[i+1][0])
    total += y
    err += abserr
print(total*8.31/200)

sig2 = s_ig2(200, 4311.57, 26.9815395, 2300)
wg = weight_gas2(sig2, 0.501, 0.821)
times2 = wg(data[:, 0]) * gas(data[:, 0])
total = 0
for i in range(len(data)-1):
    y, abserr = quad(
        senkei_hokan([data[i][0], times2[i]], [data[i+1][0], times2[i+1]]), data[i][0], data[i+1][0])
    total += y
    err += abserr
print(total*8.31/200)

pylab.plot(data[:, 0], data[:, 1])
pylab.plot(data[:, 2], data[:, 3])
pylab.plot(f, dgas)
pylab.plot(data[:, 0], dsol)
pylab.show()

