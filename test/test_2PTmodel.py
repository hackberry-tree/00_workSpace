#!/usr/bin/env python

import math

def delta(T, roh, mass, D):
    """
    融解の無次化された係数Δを算出
    """
    kb = 1.
    beta = 1/kb/T
    return 8/3*((mass*beta)**(0.5))*D*(roh**(1./3))*((6/math.pi)**(2./3))

def delta2(T, roh, D):
    """
    規格化された引数を用いた係数Δの計算式
    """
    return 8./3*((math.pi/T)**(0.5))*D*(roh**(1./3))*((6/math.pi)**(2./3))

T_star = [1.8, 1.4, 1.1, 0.9]
roh_star = [0.05, 0.40, 0.70, 0.85, 1.10]
D_star = [[0.253, 0.190, 0.153, 0.113],
          [0.203, 0.166, 0.129, 0.069],
          [0.140, 0.100, 0.080, 0.069],
          [0.087, 0.062, 0.049, 0.036],
          [0.000, 0.000, 0.000, 0.000]]

def eq34(delta, frac):
    """
    式34の左辺の値をreturn
    """
    term1 = 2*(delta**(-9./2))*(frac**(15./2))
    term2 = -6*(delta**(-3))*(frac**5)
    term3 = -1*(delta**(-3./2))*(frac**(7./2))
    term4 = 6*(delta**(-3./2))*(frac**(5./2))
    term5 = 2*frac
    term6 = -2
    return term1 + term2 + term3 + term4 + term5 + term6

fraction = [[0.936, 0.921, 0.911, 0.889],
            [0.690, 0.675, 0.647, 0.535],
            [0.534, 0.491, 0.470, 0.460],
            [0.417, 0.379, 0.358, 0.326],
            [0.0163, 0.0152, 0.0128, 0.0123]]

delta_ref = [[10.125, 8.612, 7.812, 6.399],
             [2.0224, 1.886, 1.653, 0.973],
             [0.964, 0.781, 0.703, 0.667],
             [0.529, 0.428, 0.378, 0.307],
             [1.23e-3, 1.09e-3, 8.05e-4, 7.52e-4]]

def get_f(delta, fract=0):
    """
    再帰関数からfractを求める
    """
    f = fract
    zansa = eq34(delta, fract)
    # if zansa**2 < 1e-5:
    #     return f
    # else:
    #     return get_f(delta, f + 0.001)
    while zansa < 0:
        f += 0.000001
        zansa = eq34(delta, f)
    return f, zansa

for i in range(9):
    idx = i-5
    delta = 1 * 10 ** (idx)
#    print(get_f(delta, 0))


print(get_f(0.8241459739659586))
# for i in range(4):
#     for j in range(5):
#         delta = delta2(T_star[i], roh_star[j], D_star[j][i]/roh_star[j])
#         #print(delta)
#         #print(eq34(delta_ref[j][i], fraction[j][i]))
#         print(get_f(delta_ref[j][i], 0))

