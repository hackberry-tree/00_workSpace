#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing as mp

L = 20000
proc = 2    # 8並列とする

# 各プロセスが実行する計算
def subcalc(p): # p = 0,1,...,7
    subtotal = 0

    # iの範囲を設定
    ini = L * p / proc
    fin = L * (p+1) / proc

    # 計算を実行
    for i in range(ini, fin):
        for j in range(L):
            subtotal += i * j
    return subtotal

# 8個のプロセスを用意
pool = mp.Pool(proc)

# 各プロセスに subcalc(p) を実行させる
# ここで p = 0,1,...,7
# callbackには各戻り値がlistとして格納される
callback = pool.map(subcalc, range(2))

# 各戻り値の総和を計算
total = sum(callback)

print (total)
