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

EXE_PATH = "/Users/enoki/Researches/Analysis/Codes/01_testRun/montecarlo/"

def main():
    """
    main
    """
    # a = from_file_fcc_prim_cluster(os.path.join(EXE_PATH, 'log.txt'))
    # interstitial()
    fcc_from_class()
    # tetragonal()


class CEM_Parser(object):
    """
    CEM の結果を読んで cluster, eci の情報を取り扱う
    """
    def __init__(self, clus, ecis_dict):
        if len(clus[0]) != 1:
            print("ERROR: single cluster dose not exist")

        self.clus = [clus[i] for i in sorted(ecis_dict.keys())]
        self.ecis = [ecis_dict[i] for i in sorted(ecis_dict.keys())]

    def make_pickle_quad(self, dirc):
        """
        quad site の記述に変換して pickle を出力
        """
        clus_out = [np.array([[self.conv_quad(x) for x in cl] for cl in clu])
                    for clu in self.clus[1:]]
        with open(os.path.join(dirc, "cluster.pickle"), 'wb') as wbfile:
            pickle.dump({'clusters': clus_out, 'ecis': self.ecis[1:],
                         'eci_point': self.ecis[0]}, wbfile)

    @staticmethod
    def _symm_cubic(sites):
        """
        対称操作で等価なサイトを作成する
        重複は考慮せず複製していく
        立方晶対称性の method
        n体のクラスターにつき、どこを原点にするかという並進対称操作も必要
        n体 × 48 パターンを作成する
        args:
            sites: n体クラスターの座標 [n*3 of array]
        other variables:
            nbody: 何体の cluster であるか
            sites_t: 並進対象操作した site を格納
            sites_r: 回転対象操作した site を格納
            sites_m: 反転対象操作した site を格納
            double: すべての原点位置を変えた site の
                    組み合わせを作成するため利用する一時的な変数
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

    @staticmethod
    def _symm_tetra(sites):
        """
        対称操作で等価なサイトを作成する
        重複は考慮しない
        正方晶の method
        n体のクラスターにつき、どこを原点にするかという並進対称操作も必要
        n体 × 16 パターンを作成する
        args:
            sites: n体クラスターの座標 [n*3 of array]
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

    @staticmethod
    def _rm_dup(items):
        """
        重複する要素を除外する

        以下手順
        1. 第一カラム([0, 0, 0])を除去して list 化する
        2. 各要素をソートする
        3. 重複を除外する
        4. 第一カラムに[0,0,0]を戻してreturn
        """
        # 1.
        list_i = []
        for item in items[:, 1:, :]:
            li = [list(x) for x in item]
            list_i.append(li)
        # 2.
        sorted_items = [sorted(y) for y in list_i]
        # 3.
        unique = []
        for item in sorted_items:
            if not item in unique:
                unique.append(item)
        # 4.
        out = [[[0, 0, 0]] + x for x in unique]

        return out

    @classmethod
    def from_dirc_std(cls, dirc):
        """
        dirctory を指定して
        cluster の形 (site) を log.txt から、
        ecis を ecis.txt から読む
        表記法が conventional のケース
        """
        clus = cls.parse_logtxt_std(os.path.join(dirc, "log.txt"))
        ecis = cls.parse_ecitxt(os.path.join(dirc, 'eci.txt'))
        return cls(clus, ecis)

    @classmethod
    def parse_logtxt_std(cls, path):
        """
        log.txt を parse して cluster 群を return する
        表記方法が conventional cell 用
        以下手順
        1. ファイルの読み込み
        2. meta を使って cluster の情報が記載されている位置を習得
        3. parse
        4. 対称操作で cluster を複製、重複しないもののみを抜き出す
           また原点の座標は削除して return する
        """
        # 1.
        with open(path, 'r') as rfile:
            lines = rfile.readlines()
        # 2.
        tmp_clus = []
        meta_init = re.compile(r"number of clusters in DISORDERED STATE: .*")
        init = 0
        while not meta_init.match(lines[init]):
            init += 1
        end = init + 1
        meta_end = re.compile(r"Species represented in concentration vector: .*")
        while not meta_end.match(lines[end]):
            end += 1
        # 3.
        for i in range(init+1, end-1):
            if lines[i].split()[0] == 'cluster:':
                tmp_clus.append(np.array([[float(y) for y in x.split()]
                                          for x in lines[i+1:i+4]]).T)
        # 4.
        clus = [[x[1:] for x in cls._rm_dup(cls._symm_cubic(sub_clus))]
                for sub_clus in tmp_clus]

        return clus


    @staticmethod
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


    @staticmethod
    def parse_ecitxt(fname):
        """
        null クラスターは flip に関与しないので削除する
        clus_pos と対応させるために id は -1 する
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        # num_cls = int(lines[0].split()[0])
        pos = lines.index(" 0.000000\n")
        cls_ids = []
        for i in range(1, pos):
            cls_ids += [int(x) - 1 for x in lines[i].split()]
        ecis = []
        for i in range(pos+1, len(lines)):
            ecis += [float(x) for x in lines[i].split()]
        ecis_dict = {i:x for i, x in zip(cls_ids, ecis)}
        ecis_dict.pop(-1)
        return ecis_dict


def fcc_from_class():
    """
    fcc の cluster ecis の情報を pickle にして書き出す
    return:
        dict{clusters: [(n - 2) × degen × (vertex - 1) × 4 座標]
             ecis: [(n - 2) * ecis]
             eci_point: point cluster のみの ecis}
        n-2 は null と point を除いてあることに
        vertex-1 は 原点 [0, 0, 0, 0] を除いてあることに対応する

    以下手順
    1. ecis と cluster を読み込む
    2. ecis に値があてられている cluster のみ抽出、更に pointvi  cluster を取り除く
    3.
    4. cluster のサイトの表記を quad 表記に修正
    """
    # 1.
    path = os.path.join(EXE_PATH, "AlCu/wien/TO")
    # (primitive cell 表記の場合)
    # clus = from_file_fcc_prim_cluster(os.path.join(path, "log.txt"))
    parser = CEM_Parser.from_dirc_std(path)
    parser.make_pickle_quad(path)


def fcc_from_class_orig():
    """
    fcc の cluster ecis の情報を pickle にして書き出す
    return:
        dict{clusters: [(n - 2) × degen × (vertex - 1) × 4 座標]
             ecis: [(n - 2) * ecis]
             eci_point: point cluster のみの ecis}
        n-2 は null と point を除いてあることに
        vertex-1 は 原点 [0, 0, 0, 0] を除いてあることに対応する

    以下手順
    1. ecis と cluster を読み込む
    2. ecis に値があてられている cluster のみ抽出
    3. point cluster を分離して、さらに重複する cluster を削除して
       cl, ecis_out に data を格納
    4. cluster のサイトの表記を quad 表記に修正
    """
    # 1.
    path = os.path.join(EXE_PATH, "AlCu/wien/TO")
    # (primitive cell 表記の場合)
    # clus = from_file_fcc_prim_cluster(os.path.join(path, "log.txt"))
    parser = CEM_Parser.from_dirc_std(path)
    clus = parser.out
    ecis = parser.ecis
    # 2.
    tmp = [clus[i] for i in ecis[0]]
    # 3.
    ecis_out = []
    cl = []
    if len(clus[0]) != 1:
        print("ERROR: single cluster dose not exist")
    for i in range(len(ecis[0]) - 1):
        uni = [x[1:] for x in rm_dup(symm_cubic(tmp[i+1]))]
        ecis_out.append(ecis[1][ecis[0][i+1]])
        cl.append(uni)
    # 4.
    cl_out = [np.array([[conv_quad(x) for x in cl2] for cl2 in cl1])
              for cl1 in cl]
    with open(os.path.join(path, "cluster.pickle"), 'wb') as wbfile:
        pickle.dump({'clusters': cl_out, 'ecis': ecis_out,
                     'eci_point': ecis[1][0]}, wbfile)
    return cl_out


def symm_cubic(sites):
    """
    対称操作で等価なサイトを作成する
    重複は考慮しない
    立方晶の method
    n体のクラスターにつき、どこを原点にするかという並進対称操作も必要
    n体 × 48 パターンを作成する
    args:
        sites: n体クラスターの座標 [n*3 of array]
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
        sites: n体クラスターの座標 [n*3 of array]
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
    重複する要素を除外する

    以下手順
    1. 第一カラム([0, 0, 0])を除去して list 化する
    2. 各要素をソートする
    3. 重複を除外する
    4. 第一カラムに[0,0,0]を戻してreturn
    """
    # 1.
    list_i = []
    for item in items[:, 1:, :]:
        li = [list(x) for x in item]
        list_i.append(li)
    # 2.
    sorted_items = [sorted(y) for y in list_i]
    # 3.
    unique = []
    for item in sorted_items:
        if not item in unique:
            unique.append(item)
    # 4.
    out = [[[0, 0, 0]] + x for x in unique]

    return out


def from_file_fcc_cluster(fname):
    """
    cluster の形 (site) を log.txt から読み込む
    """
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


def from_file_fcc_prim_cluster(fname):
    """
    cluster の site を log.txt から読み込む
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


def from_file_int_cluster(fname):
    """
    cluster の site を log.txt から読み込む
    i-s 用
    """
    with open(fname, 'r') as rfile:
        lines = rfile.readlines()
    out = []
    for i in range(len(lines)):
        if lines[i].split()[0] == 'cluster:':
            out.append(np.array([[float(y) for y in x.split()]
                                 for x in lines[i+1:i+4]]).T)
    return out


def from_file_int2(fname):
    """
    cluster のサイトの置換元素が A/C どちらかを読み取る
    """
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
    # num_cls = int(lines[0].split()[0])
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
    """
    対象操作を行った後、重複を削除
    """
    out = []
    for i in range(len(clus[0])):
        dups = symm_cube_int(clus[0][i], clus[1][i])
        uni = rm_dup(dups)
        out.append(
            np.array([[conv_octa(site) for site in cluster] for cluster in uni]))
    return out


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
    """
    i-s 用の構造読み込み
    """
    clus = from_file_int_cluster(os.path.join(EXE_PATH, 'log_int_Ti.txt'))
    clus_pos = from_file_int2(os.path.join(EXE_PATH, 'log_int_Ti.txt'))
    ecis = ecis_from_file(os.path.join(EXE_PATH, 'eci_int_Ti_onlyNN.txt'))
    tmp = ext_used(clus, clus_pos, ecis)
    a = get_all(tmp['a'])
    c = get_all(tmp['c'])
    eci_a = tmp['a'][2]
    eci_c = tmp['c'][2]
    with open(os.path.join(EXE_PATH, "cluster_int_Ti_onlyNN.pickle"), 'wb') as wbfile:
        pickle.dump({'sub_clusters': c, 'int_clusters': a,
                     'sub_ecis': eci_c, 'int_ecis': eci_a}, wbfile)


def fcc():
    """
    fcc の cluster ecis の情報を pickle にして書き出す
    return:
        dict{clusters: [(n - 2) × degen × (vertex - 1) × 4 座標]
             ecis: [(n - 2) * ecis]
             eci_point: point cluster のみの ecis}
        n-2 は null と point を除いてあることに
        vertex-1 は 原点 [0, 0, 0, 0] を除いてあることに対応する

    以下手順
    1. ecis と cluster を読み込む
    2. ecis に値があてられている cluster のみ抽出
    3. point cluster を分離して、さらに重複する cluster を削除して
       cl, ecis_out に data を格納
    4. cluster のサイトの表記を quad 表記に修正
    """
    # 1.
    path = os.path.join(EXE_PATH, "AlCu/wien/TO")
    # (primitive cell 表記の場合)
    # clus = from_file_fcc_prim_cluster(os.path.join(path, "log.txt"))
    clus = from_file_fcc_cluster(os.path.join(path, "log.txt"))
    ecis = ecis_from_file(os.path.join(path, "eci.txt"))
    # 2.
    tmp = [clus[i] for i in ecis[0]]
    # 3.
    ecis_out = []
    cl = []
    if len(clus[0]) != 1:
        print("ERROR: single cluster dose not exist")
    for i in range(len(ecis[0]) - 1):
        uni = [x[1:] for x in rm_dup(symm_cubic(tmp[i+1]))]
        ecis_out.append(ecis[1][ecis[0][i+1]])
        cl.append(uni)
    # 4.
    cl_out = [np.array([[conv_quad(x) for x in cl2] for cl2 in cl1])
              for cl1 in cl]
    with open(os.path.join(path, "cluster.pickle"), 'wb') as wbfile:
        pickle.dump({'clusters': cl_out, 'ecis': ecis_out,
                     'eci_point': ecis[1][0]}, wbfile)
    print(path)
    return cl_out


def tetragonal():
    """
    dict{clusters: [(n - 2) × degen × (vertex - 1) × 4 座標]
         ecis: [(n - 2) * ecis]
         eci_point: point cluster のみの ecis}
    n-2 は null と point を除いてあることに
    vertex-1 は 原点 [0, 0, 0, 0] を除いてあることに対応する
    """
    path = os.path.join(EXE_PATH, "AlCu/voldep/4.0/")
    # point 以上のクラスター修得
    clus = from_file_fcc_prim_cluster(os.path.join(path, "log.txt"))
    cl = []
    # point 以上の eci 修得
    ecis = ecis_from_file(os.path.join(path, 'eci.txt'))
    tmp = [clus[i] for i in ecis[0]]
    ecis_out = []
    cl = []
    if len(clus[0]) != 1:
        print("ERROR: single cluster dose not exist")
    for i in range(len(ecis[0]) - 1):
        uni = [x[1:] for x in rm_dup(symm_tetra(tmp[i+1]))]
        ecis_out.append(ecis[1][ecis[0][i+1]])
        cl.append(uni)
    conv_cl = [np.array([[conv_quad(x) for x in cl2] for cl2 in cl1])
               for cl1 in cl]
    cl = conv_cl
    with open(os.path.join(path, "cluster.pickle"), 'wb') as wbfile:
        pickle.dump({'clusters': cl, 'ecis': ecis_out,
                     'eci_point': ecis[1][0]}, wbfile)


if __name__ == '__main__':
    main()

