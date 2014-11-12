#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
x yデータを線で結ぶためにソートすることを考える
円のようにxやyが+から-へ戻るような場合でもうまく処理したいスクリプト
1: xで昇順ソートする
2: 指定のdy以内のデータ収集
3: 最も近い距離のものを選択
"""

import copy
import numpy as np
import pylab

__date__ = "Sep 4 2014"

PATH = "/Users/enoki/Desktop/cvm02/re_plot_cvm/dat_00.txt"

def main():
    """main"""
    with open(PATH, 'r') as rfile:
        lines = rfile.readlines()
    data = [[float(x.split()[1]), float(x.split()[0])] for x in lines]
    data.sort(key=lambda x: x[0])
    data = set_for_list(data)
    data.sort(key=lambda x: x[1])

    tail = data.pop(0)
    within = search_within_dy(2, tail, data)
    out = [tail]
    while within:
        tail = search_nearlest(tail, within, 800)
        out.append(tail)
        data.remove(tail)
        within = search_within_dy(50, tail, data)
    for pt in out:
        print("{0[0]} {0[1]}".format([str(pt[0]), str(pt[1])]))

    datax = [x[0] for x in out]
    datay = [x[1] for x in out]


    pylab.plot(datax, datay)


    data = [[float(x.split()[3]), float(x.split()[2])] for x in lines]
    data.sort(key=lambda x: x[0])
    data = set_for_list(data)
    data.sort(key=lambda x: x[1])

    tail = data.pop(0)
    within = search_within_dy(2, tail, data)
    out = [tail]
    while within:
        tail = search_nearlest(tail, within, 800)
        out.append(tail)
        data.remove(tail)
        within = search_within_dy(1000, tail, data)
    for pt in out:
        print("{0[0]} {0[1]}".format([str(pt[0]), str(pt[1])]))

    datax = [x[0] for x in out]
    datay = [x[1] for x in out]


    pylab.plot(datax, datay)
    pylab.show()

def collect_path(original, tail, yovx=800):
    pass

def search_within_dy(val, p0, ps):
    return ps
    out = []
    for pt in ps:
        if (p0[1] - pt[1]) ** 2 <= val ** 2:
            out.append(pt)
    return out


def search_nearlest(p0, ps, yovx):
    tmp = [[get_distance(p0, x, yovx), x] for x in ps]
    tmp.sort(key=lambda x: x[0])
    return tmp[0][1]

def get_distance(p1, p2, yovx):
    """
    yova: xとyの比率を入れる グラフのy/x
    """
    dp = (np.array(p1) - np.array(p2))
    dp[0] *= yovx
    return np.linalg.norm(dp)

def set_for_list(in_list):
    """
    リストのリストで重複要素を削除
    データは予めソートしておく必要がある
    """
    # 重複要素削除 リストデータなのでsetが使えない
    i = 0
    out_list = copy.deepcopy(in_list)
    while len(out_list) > i + 1:
        if out_list[i] == out_list[i+1]:
            out_list.pop(i)
        else:
            i += 1
    return out_list



if __name__ == "__main__":
    main()
