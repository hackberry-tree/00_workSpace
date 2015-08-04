#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
from itertools import combinations, chain, product
import numpy as np
import os
import pickle
import math
import re

# TEST_PATH = "/Users/enoki"
TEST_PATH = "/Users/enoki/Researches/Analysis/Codes/01_testRun/montecarlo/"

def symm_cube(sites):
    """
    対称操作で等価なサイトを作成する
    重複は考慮しない
    立方晶の method
    n体のクラスターにつき、どこを原点にするかという並進対称操作も必要
    n体 × 48 パターンを作成する
    args:
        sites: n体クラスターの座標[n*3 of array]
    """
    nbody = sites.shape[0]
    sites_t = []
    double = np.r_[sites, sites]
    # 先頭のサイトは常に[0, 0, 0]
    for i in range(len(sites)):
        sites_t.append((double - sites[i])[i:i+nbody])
    sites_t = np.array(sites_t)

    symm_r = np.array([[0, 1, 2], [1, 0, 2], [0, 2, 1],
                       [1, 2, 0], [2, 0, 1], [2, 1, 0]])
    sites_r = []
    for site in sites_t:
        for symm in symm_r:
            sites_r.append(site[:, symm])
    sites_r = np.array(sites_r)

    symm_m = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1],
                       [-1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])
    return (sites_r * symm_m.reshape(8, 1, 1, 3)).reshape(48*nbody, nbody, 3)


def symm_tetra(sites):
    """
    対称操作で等価なサイトを作成する
    重複は考慮しない
    正方晶の method
    n体のクラスターにつき、どこを原点にするかという並進対称操作も必要
    n体 × 16 パターンを作成する
    args:
        sites: n体クラスターの座標[n*3 of array]
    """
    nbody = sites.shape[0]
    sites_t = []
    double = np.r_[sites, sites]
    # 先頭のサイトは常に[0, 0, 0]
    for i in range(len(sites)):
        sites_t.append((double - sites[i])[i:i+nbody])
    sites_t = np.array(sites_t)

    symm_r = np.array([[0, 1, 2], [1, 0, 2]])
    sites_r = []
    for site in sites_t:
        for symm in symm_r:
            sites_r.append(site[:, symm])
    sites_r = np.array(sites_r)

    symm_m = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1],
                       [-1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])
    return (sites_r * symm_m.reshape(8, 1, 1, 3)).reshape(16*nbody, nbody, 3)


def rm_dup(items):
    """
    0. リスト化する
    0.1. 第一カラム([0, 0, 0])を除去
    1. 各要素をソートする
    2. 重複を除外する
    3. [0,0,0]を戻してreturn
    """
    list_i = []
    for item in items[:,1:,:]:
        li = [list(x) for x in item]
        list_i.append(li)
    sorted_items = [sorted(y) for y in list_i]
    unique = []
    for item in sorted_items:
        if not item in unique:
            unique.append(item)
    out = [[[0, 0, 0]] + x for x in unique]
    return out


def from_file(fname):
    with open(fname, 'r') as rfile:
        lines = rfile.readlines()
    out = []

    meta_init = re.compile(r"number of clusters in DISORDERED STATE: .*")
    init = 0
    while not meta_init.match(lines[init]):
        init += 1
    end = init + 1
    meta_end = re.compile(r"Species represented in concentration vector: .*")
    while not meta_end.match(lines[end]):
        end += 1

    for i in range(init+1, end-1):
        if lines[i].split()[0] == 'cluster:':
            out.append(np.array([[float(y) for y in x.split()]
                       for x in lines[i+1:i+4]]).T)
    return out


def from_file2(fname):
    """
    表記法が fcc primitive のケース
    """
    with open(fname, 'r') as rfile:
        lines = rfile.readlines()
    out = []
    a = np.array([[1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2]])

    meta_init = re.compile(r"number of clusters in DISORDERED STATE: .*")
    init = 0
    while not meta_init.match(lines[init]):
        init += 1
    end = init + 1
    meta_end = re.compile(r"Species represented in concentration vector: .*")
    while not meta_end.match(lines[end]):
        end += 1

    for i in range(init+1, end-1):
        if lines[i].split()[0] == 'cluster:':
            out.append(np.array([[float(y) for y in x.split()]
                       for x in lines[i+1:i+4]]).T)
    out2 = []
    for clus in out:
        tmpc = []
        for s in clus:
            new_s = (a * s.T).sum(-1).T
            tmpc.append(new_s)
        out2.append(np.array(tmpc))
    print(out)
    print(out2)
    # for i in range(len(lines)):
    #     if lines[i].split()[0] == 'cluster:':
    #         out.append(np.array([list((a*np.array([float(y) for y in x.split()]).T).sum(-1).T)
    #                    for x in lines[i+1:i+4]]).T)
    return out2


def conv_quad(site):
    """
    QuadSiteの表記に変換する
    """
    integ = []
    decim = []
    for coord in site:
         integ.append(math.floor(coord))
         decim.append(coord - math.floor(coord))
    site_id = {(0., 0., 0.): 0, (0., 0.5, 0.5): 1,
               (0.5, 0., 0.5): 2, (0.5, 0.5, 0.): 3}[tuple(decim)]
    integ.append(site_id)
    return integ


def from_file_int(fname):
    with open(fname, 'r') as rfile:
        lines = rfile.readlines()
    out = []
    for i in range(len(lines)):
        if lines[i].split()[0] == 'cluster:':
            out.append(np.array([[float(y) for y in x.split()]
                       for x in lines[i+1:i+4]]).T)
    return out


def from_file_int2(fname):
    with open(fname, 'r') as rfile:
        lines = rfile.readlines()
    out = []
    i = 0
    while lines[i] != "      # disord  ordrd          iants\n":
        i += 1
    for line in lines[i+1:]:
        out.append(line.split()[6])
    pos_c = [[i for i in range(len(x)) if x[i] == 'C'] for x in out]
    in_c = [i for i in range(len(pos_c)) if pos_c[i]]
    pos_a = [[i for i in range(len(x)) if x[i] == 'A'] for x in out]
    in_a = [i for i in range(len(pos_a)) if pos_a[i]]
    return in_c, pos_c, in_a, pos_a


def ecis_from_file(fname):
    """
    null クラスターは flip に関与しないので削除する
    clus_pos と対応させるために id は -1 する
    """
    with open(fname, 'r') as rfile:
        lines = rfile.readlines()
    num_cls = int(lines[0].split()[0])
    pos = lines.index(" 0.000000\n")
    cls_ids = []
    for i in range(1, pos):
        cls_ids += [int(x) - 1 for x in lines[i].split()]
    ecis = []
    for i in range(pos+1, len(lines)):
        ecis += [float(x) for x in lines[i].split()]
    clus = {i:x for i, x in zip(cls_ids, ecis)}
    clus.pop(-1)
    cls_ids.sort()
    ecis = [clus[i] for i in cls_ids[1:]]
    return cls_ids[1:], clus


def ext_used(cluster, clpos, ecis):
    """
    eci_ids に表記のある cluster のみを抜き出す
    """
    c_id = [i for i in ecis[0] if i in clpos[0]]
    a_id = [i for i in ecis[0] if i in clpos[2]]
    c_cls = [cluster[i] for i in c_id]
    c_ord = [clpos[1][i] for i in c_id]
    c_ecis = [ecis[1][i] for i in c_id]
    a_cls = [cluster[i] for i in a_id]
    a_ord = [clpos[3][i] for i in a_id]
    a_ecis = [ecis[1][i] for i in a_id]
    return {'c':[c_cls, c_ord, c_ecis], 'a':[a_cls, a_ord, a_ecis]}


def ext_used_single(cluster, ecis):
    """
    eci_ids に表記のある cluster のみを抜き出す
    """
    c_cls = [cluster[i] for i in ecis[0]]
    return [c_cls, ecis]


def symm_cube_int(sites, pos):
    """
    対称操作で等価なサイトを作成する
    重複は考慮しない
    立方晶の method
    n体中の等価なサイトmにつき、どこを原点にするかという並進対称操作も必要
    m × 48 パターンを作成する
    args:
        sites: n体クラスターの座標 [n*3 of array]
        pos: 原点と等価なサイト位置 [m of list]
    """
    nbody = sites.shape[0]
    sites_t = []
    double = np.r_[sites, sites]
    # 先頭のサイトは常に[0, 0, 0]
    for i in pos:
        sites_t.append((double - sites[i])[i:i+nbody])
    sites_t = np.array(sites_t)

    symm_r = np.array([[0, 1, 2], [1, 0, 2], [0, 2, 1],
                       [1, 2, 0], [2, 0, 1], [2, 1, 0]])
    sites_r = []
    for site in sites_t:
        for symm in symm_r:
            sites_r.append(site[:, symm])
    sites_r = np.array(sites_r)

    symm_m = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1],
                       [-1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])

    return (sites_r * symm_m.reshape(8, 1, 1, 3)).reshape(48*len(pos), nbody, 3)


def get_all(clus):
    out = []
    for i in range(len(clus[0])):
        dups = symm_cube_int(clus[0][i], clus[1][i])
        uni = rm_dup(dups)
        out.append(
            np.array([[conv_octa(site) for site in cluster] for cluster in uni]))
    return(out)


def conv_octa(site):
    """
    OctaSiteの表記に変換する
    """
    integ = []
    decim = []
    for coord in site:
         integ.append(math.floor(coord))
         decim.append(coord - math.floor(coord))
    site_id = {(0.0, 0.0, 0.0): 0, (0.0, 0.5, 0.5): 1,
               (0.5, 0.0, 0.5): 2, (0.5, 0.5, 0.0): 3,
               (0.5, 0.0, 0.0): 4, (0.0, 0.5, 0.0): 5,
               (0.0, 0.0, 0.5): 6, (0.5, 0.5, 0.5): 7}[tuple(decim)]
    integ.append(site_id)
    return integ


def interstitial():
    clus = from_file_int(os.path.join(TEST_PATH, 'log_int_Ti.txt'))
    clus_pos = from_file_int2(os.path.join(TEST_PATH, 'log_int_Ti.txt'))
    ecis = ecis_from_file(os.path.join(TEST_PATH, 'eci_int_Ti_onlyNN.txt'))
    tmp = ext_used(clus, clus_pos, ecis)
    a = get_all(tmp['a'])
    c = get_all(tmp['c'])
    eci_a = tmp['a'][2]
    eci_c = tmp['c'][2]
    with open(os.path.join(TEST_PATH, "cluster_int_Ti_onlyNN.pickle"), 'wb') as wbfile:
        pickle.dump({'sub_clusters': c, 'int_clusters': a,
                     'sub_ecis': eci_c, 'int_ecis': eci_a}, wbfile)


def fcc():
    """
    dict{clusters: [(n - 2) × degen × (vertex - 1) × 4 座標]
         ecis: [(n - 2) * ecis]
         eci_point: point cluster のみの ecis}
    n-2 は null と point を除いてあることに
    vertex-1 は 原点 [0, 0, 0, 0] を除いてあることに対応する
    """
    path = os.path.join(TEST_PATH, "AlCu/wien/
    # point 以上のクラスター修得
    print(path)
    # clus = from_file2(os.path.join(path, "log.txt"))
    clus = from_file(os.path.join(path, "log.txt"))
    cl = []
    # point 以上の eci 修得
    ecis = ecis_from_file(os.path.join(path, 'eci.txt'))
    tmp = ext_used_single(clus, ecis)
    ecis_out = []
    cl = []
    if len(clus[0]) != 1:
        print("ERROR: single cluster dose not exist")
    for i in range(len(ecis[0]) - 1):
        uni = [x[1:] for x in rm_dup(symm_cube(tmp[0][i+1]))]
        ecis_out.append(ecis[1][ecis[0][i+1]])
        cl.append(uni)
    conv_cl = [np.array([[conv_quad(x) for x in cl2] for cl2 in cl1])
               for cl1 in cl]
    cl = conv_cl
    with open(os.path.join(path, "cluster.pickle"), 'wb') as wbfile:
        pickle.dump({'clusters': cl, 'ecis': ecis_out,
                     'eci_point': ecis[1][0]}, wbfile)
    print(path)


def tetragonal():
    """
    dict{clusters: [(n - 2) × degen × (vertex - 1) × 4 座標]
         ecis: [(n - 2) * ecis]
         eci_point: point cluster のみの ecis}
    n-2 は null と point を除いてあることに
    vertex-1 は 原点 [0, 0, 0, 0] を除いてあることに対応する
    """
    path = os.path.join(TEST_PATH, "AlCu/voldep/4.0/")
    # point 以上のクラスター修得
    clus = from_file2(os.path.join(path, "log.txt"))
    cl = []
    # point 以上の eci 修得
    ecis = ecis_from_file(os.path.join(path, 'eci.txt'))
    tmp = ext_used_single(clus, ecis)
    ecis_out = []
    cl = []
    if len(clus[0]) != 1:
        print("ERROR: single cluster dose not exist")
    for i in range(len(ecis[0]) - 1):
        uni = [x[1:] for x in rm_dup(symm_tetra(tmp[0][i+1]))]
        ecis_out.append(ecis[1][ecis[0][i+1]])
        cl.append(uni)
    conv_cl = [np.array([[conv_quad(x) for x in cl2] for cl2 in cl1])
               for cl1 in cl]
    cl = conv_cl
    with open(os.path.join(path, "cluster.pickle"), 'wb') as wbfile:
        pickle.dump({'clusters': cl, 'ecis': ecis_out,
                     'eci_point': ecis[1][0]}, wbfile)

# a = from_file2(os.path.join(TEST_PATH, 'log.txt'))
# interstitial()
fcc()
# tetragonal()
