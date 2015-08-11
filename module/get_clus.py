#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
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
    bcci()
    fcci()
    fcc()
    fcc_prim()


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

    @classmethod
    def get_unique(cls, clus):
        """
        対象操作を行って
        重複を削除
        サイトの表記を NaCl に変換
        """
        unique_clus = []
        for i in range(len(clus[0])):
            sym_dups = cls._symm_cubic(clus[0][i], clus[1][i])
            sym_uni = cls._rm_dup(sym_dups)
            unique_clus.append(
                np.array([[cls.conv_NaCl(site) for site in cluster]
                          for cluster in sym_uni]))
        return unique_clus


    @staticmethod
    def _symm_cubic(sites, equiv_pos=None):
        """
        対称操作で等価なサイトを作成する
        重複は考慮しない
        立方晶の method
        n体中の等価なサイトmにつき、どこを原点にするかという並進対称操作も必要
        m × 48 パターンを作成する
        副格子系は equiv_pos が必要
        args:
            sites: n体クラスターの座標 [n*3 of array]
            equiv_pos: 原点と等価なサイト位置 [m of list]
        """
        nbody = len(sites)
        sites_t = []
        double = np.r_[sites, sites]
        # 先頭のサイトは常に[0, 0, 0]
        if not equiv_pos:
            equiv_pos = range(nbody)
        for i in equiv_pos:
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
        return (sites_r * symm_m.reshape(8, 1, 1, 3)).reshape(48*len(equiv_pos),
                                                              nbody, 3)

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
    def from_dirc_prim_fcc(cls, dirc):
        """
        dirctory を指定して
        cluster の形 (site) を log.txt から、
        ecis を ecis.txt から読む
        表記法が conventional のケース
        """
        clus = cls.parse_logtxt_prim_fcc(os.path.join(dirc, "log.txt"))
        ecis = cls.parse_ecitxt(os.path.join(dirc, 'eci.txt'))
        return cls(clus, ecis)

    @classmethod
    def from_dirc_2sub(cls, dirc):
        """
        dirctory を指定して
        2副格子のモデルを読み込む
        """
        clus, pos = cls.parse_logtxt_2sub(os.path.join(dirc, "log.txt"))
        ecis = cls.parse_ecitxt(os.path.join(dirc, 'eci.txt'))
        clus_dict = cls.extract_used_cluster(clus, pos, ecis)
        a = cls.get_unique(clus_dict['a'])
        c = cls.get_unique(clus_dict['c'])
        eci_a = clus_dict['a'][2]
        eci_c = clus_dict['c'][2]
        with open(os.path.join(dirc, "cluster.pickle"), 'wb') as wbfile:
            pickle.dump({'sub_clusters': c, 'int_clusters': a,
                         'sub_ecis': eci_c, 'int_ecis': eci_a}, wbfile)


    @staticmethod
    def extract_used_cluster(cluster, clpos, ecis):
        """
        eci_ids に表記のある cluster のみを抜き出す
        """
        c_id = [i for i in sorted(ecis.keys()) if i in clpos[0]]
        a_id = [i for i in sorted(ecis.keys()) if i in clpos[2]]
        try:
            c_cls = [cluster[i] for i in c_id]
        except IndexError:
            print(len(cluster))
            print(c_id)
        c_ord = [clpos[1][i] for i in c_id]
        c_ecis = [ecis[i] for i in c_id]
        a_cls = [cluster[i] for i in a_id]
        a_ord = [clpos[3][i] for i in a_id]
        a_ecis = [ecis[i] for i in a_id]
        return {'c':[c_cls, c_ord, c_ecis], 'a':[a_cls, a_ord, a_ecis]}



    @staticmethod
    def parse_logtxt_2sub(fname):
        """
        cluster の site を log.txt から読み込む
        i-s 用
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()

        meta_init = re.compile(r"number of clusters in DISORDERED STATE: .*")
        init = -1
        while not meta_init.match(lines[init]):
            init -= 1
        midle = init + 1
        meta_midle = re.compile(r"CF:     clust# clust#   dcrtn# #var- degen  dcrtn.*")
        while not meta_midle.match(lines[midle]):
            midle += 1
        end = midle + 1
        meta_end = re.compile(r"\*ECLI\*.*")
        while not meta_end.match(lines[end]):
            end += 1

        clus = []
        for i in range(init+1, midle-1):
            if lines[i].split()[0] == 'cluster:':
                clus.append(np.array([[float(y) for y in x.split()]
                                      for x in lines[i+1:i+4]]).T)
        # print(len(clus))

        tmp = [line.split()[6] for line in lines[midle+2: end-4]]
        # print(tmp[0], tmp[-1]) # 最初と最後の pos が一致するか check
        pos_c = [[i for i in range(len(x)) if x[i] == 'C'] for x in tmp]
        in_c = [i for i in range(len(tmp)) if pos_c[i]]
        pos_a = [[i for i in range(len(x)) if x[i] == 'A'] for x in tmp]
        in_a = [i for i in range(len(tmp)) if pos_a[i]]

        return clus, (in_c, pos_c, in_a, pos_a)

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
        meta_init = re.compile(r"number of clusters in DISORDERED STATE: .*")
        init = 0
        while not meta_init.match(lines[init]):
            init += 1
        end = init + 1
        meta_end = re.compile(r"Species represented in concentration vector: .*")
        while not meta_end.match(lines[end]):
            end += 1
        # 3.
        tmp_clus = []
        for i in range(init+1, end-1):
            if lines[i].split()[0] == 'cluster:':
                tmp_clus.append(np.array([[float(y) for y in x.split()]
                                          for x in lines[i+1:i+4]]).T)
        # 4.
        clus = [[x[1:] for x in cls._rm_dup(cls._symm_cubic(sub_clus))]
                for sub_clus in tmp_clus]

        return clus

    @classmethod
    def parse_logtxt_prim_fcc(cls, path):
        """
        log.txt を parse して cluster 群を return する
        表記方法が primitive fcc cell 用
        以下手順
        1. ファイルの読み込み
        2. meta を使って cluster の情報が記載されている位置を習得
        3. parse
        4. デカルト座標に変換
        5. 対称操作で cluster を複製、重複しないもののみを抜き出す
           また原点の座標は削除して return する
        """
        # 1.
        with open(path, 'r') as rfile:
            lines = rfile.readlines()
        # 2.
        meta_init = re.compile(r"number of clusters in DISORDERED STATE: .*")
        init = 0
        while not meta_init.match(lines[init]):
            init += 1
        end = init + 1
        meta_end = re.compile(r"Species represented in concentration vector: .*")
        while not meta_end.match(lines[end]):
            end += 1
        # 3.
        prim = []
        for i in range(init+1, end-1):
            if lines[i].split()[0] == 'cluster:':
                prim.append(np.array([[float(y) for y in x.split()]
                                      for x in lines[i+1:i+4]]).T)
        axis = np.array([[1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2]])
        # 4.
        tmp_clus = [np.array([(axis * cood.T).sum(-1).T for cood in clus])
                    for clus in prim]
        print(prim)
        print(tmp_clus)

        # 5.
        clus = [[x[1:] for x in cls._rm_dup(cls._symm_cubic(sub_clus))]
                for sub_clus in tmp_clus]

        return clus

    @staticmethod
    def conv_quad(site):
        """
        QuadSiteの表記に変換する
        dcimal to integer
        """
        return
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
    def conv_NaCl(site):
        """
        NaCl site の表記に変換する
        dcimal to integer
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

    @classmethod
    def exe_fcc(cls, path):
        """
        fcc の cluster, ecis の情報を pickle にして書き出す
        """
        cls.from_dirc_std(path).make_pickle_quad(path)

    @classmethod
    def exe_fcc_prim(cls, path):
        """
        fcc_prim の cluster, ecis の情報を pickle にして書き出す
        """
        cls.from_dirc_prim_fcc(path).make_pickle_quad(path)

    @classmethod
    def exe_fcci(cls, path):
        """
        fcci の cluster, ecis の情報を pickle にして書き出す
        """
        cls.from_dirc_2sub(path)

    @classmethod
    def exe_bcci(cls, path):
        """
        bcci の cluster, ecis の情報を pickle にして書き出す
        """
        cls.from_dirc_2sub(path)




def fcc():
    """
    fcc の cluster ecis の情報を pickle にして書き出す
    """
    path = os.path.join(EXE_PATH, "AlCu/wien/TO")
    CEM_Parser.exe_fcc(path)

def fcc_prim():
    """
    fcc の cluster ecis の情報を pickle にして書き出す
    """
    path = os.path.join(EXE_PATH, "AlCu/tetra_2R2N/")
    parser = CEM_Parser.from_dirc_prim_fcc(path)
    parser.make_pickle_quad(path)

def fcci():
    """
    fcci 用の構造読み込み
    """
    path = os.path.join(EXE_PATH, 'i-s', 'fcci', 'TiC')
    CEM_Parser.from_dirc_2sub(path)

def bcci():
    """
    bcci 用の eci, cluster 読み込み
    """
    path = os.path.join(EXE_PATH, 'i-s', 'bcci', 'CrC')
    CEM_Parser.from_dirc_2sub(path)


if __name__ == '__main__':
    main()

