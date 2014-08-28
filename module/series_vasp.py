#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
VASPの一連の計算ファイルを準備するモジュール
INCARのパラメータはPOSCARに依存するので、
シリーズを作成する際は、常にPOSCARから先に変更する必要がある
"""
import os
import copy
import vaspy
from commopy import Bash, Cabinet


class Series(object): #pylint: disable=R0903
    """
    Series Class
    ToDo: kp_rxとkp_socのパラメータに初期値を用意してしまうと
          vaspyでの変更が反映されないので後々困りそう...要修正
    """
    def __init__(self, path_poscar, out_path):
        poscar = vaspy.Poscar(path_poscar)
        incar = vaspy.IncarReadPoscar(path_poscar)
        self.series = [{'poscar': poscar, 'path': out_path,
                        'incar': incar, 'kp_rx': 0.15, 'kp_soc': 0.11}]


class PoscarMixin(object):
    """
    Alt parameters of Poscar file.
    Control below attribute.
    """
    series = []

    def __set_param_poscar(self, method, key, param_list):
        """
        vaspy.Poscarのmethodを使ってパラメータの変更
        pathの追加
        parameterの追加
        またPOSCARに依って依存するincarのパラメータの更新を行う
        """
        out_list = []
        for posc_base in self.series:
            for param in param_list:
                posc_new = copy.deepcopy(posc_base)
                getattr(posc_new['poscar'], method)(param)
                posc_new.update({key: param})
                dst_path = os.path.join(posc_new['path'],
                                        "{0}_{1}".format(key, param))
                incar = vaspy.IncarLoadPoscarObj(posc_new['poscar'])
                posc_new.update({'path': dst_path, 'incar': incar})
                out_list.append(posc_new)
        self.series = out_list

    def set_cova(self, cova_list):
        """
        cova_listに基づいてc/aを変更
        """
        self.__set_param_poscar('alt_c_over_a', 'cova', cova_list)

    def set_scale(self, scale_list):
        """
        lattice constant (scale)を変える
        """
        self.__set_param_poscar('alt_cell_scale', 'scale', scale_list)

    def set_volume(self, volume_list):
        """
        cell volume(scale)を変える
        """
        self.__set_param_poscar('alt_cell_volume', 'volume', volume_list)

    def set_elements(self, elements_list):
        """
        elementsを変更する
        num_atomsは一定
        """
        out_list = []
        for posc_base in self.series:
            for elements in elements_list:
                posc_new = copy.deepcopy(posc_base)
                posc_new['poscar'].elements = elements
                posc_new.update({'elem': elements})
                incar = vaspy.IncarLoadPoscarObj(posc_new['poscar'])
                formula = incar.get_formula()
                dst_path = os.path.join(posc_new['path'],
                                        "{0}_{1}".format('elem', formula))
                posc_new.update({'path': dst_path, 'incar': incar})
                out_list.append(posc_new)
        self.series = out_list


class IncarMixin(object):
    """
    Alt parameters of INCAR file.
    Control below attribute.
    KPOINTSもここで変更する
    """
    series = []

    def set_incar_tag(self, key, param_list):
        """
        incarパラメータの変更
        """
        out_list = []
        for inc_base in self.series:
            for param in param_list:
                inc_new = copy.deepcopy(inc_base)
                inc_new['incar'][key] = param
                inc_new.update({key: param})
                dst_path = os.path.join(inc_new['path'],
                                        "{0}_{1}".format(key, param))
                inc_new.update({'path': dst_path})
                out_list.append(inc_new)
        self.series = out_list

    def set_incar_fixedtag(self, key, param_list):
        """
        incarのfixedtagを変更
        すべてのINCARファイルで共通になる
        """
        out_list = []
        for inc_base in self.series:
            for param in param_list:
                inc_new = copy.deepcopy(inc_base)
                inc_new['incar'].fixed_tag.update({key: param})
                inc_new.update({key: param})
                dst_path = os.path.join(inc_new['path'],
                                        "{0}_{1}".format(key, param))
                inc_new.update({'path': dst_path})
                out_list.append(inc_new)
        self.series = out_list


class KpointsMixin(object):
    """
    Alt Kpoints
    """
    series = []

    def set_kp_relax(self, param_list):
        """
        KPOINTSのdQの変更
        pathの追加
        parameterの追加
        またPOSCARに依って依存するincarのパラメータの更新を行う
        """
        out_list = []
        for posc_base in self.series:
            for param in param_list:
                posc_new = copy.deepcopy(posc_base)
                posc_new.update({'kp_rx': param})
                dst_path = os.path.join(posc_new['path'],
                                        "{0}_{1}".format('kp_rx', param))
                incar = vaspy.IncarLoadPoscarObj(posc_new['poscar'])
                posc_new.update({'path': dst_path, 'incar': incar})
                out_list.append(posc_new)
        self.series = out_list

    def set_kp_soc(self, param_list):
        """
        socのk点を指定
        """
        pass


class Produce(Series, PoscarMixin, IncarMixin, KpointsMixin):
    """
    seriesを作成する
    """
    def __init__(self, path_poscar, out_path):
        Series.__init__(self, path_poscar, out_path)

    def make_files(self):
        """
        ディレクトリを作成して、ファイルを作成
        """
        self.mkdirs()
        self.__write_poscars()
        self.__make_inputs()

    def __write_poscars(self):
        """
        self.seriesのPOSCARをpath中に展開
        """
        for param in self.series:
            path = param['path']
            dst_path = os.path.join(path, 'POSCAR')
            param['poscar'].write_poscar(dst_path)

    def __make_inputs(self):
        """
        self.seriesのINCARをpath中に展開
        """
        for param in self.series:
            dst_path = param['path']
            incar = param['incar']
            kp_rx = param['kp_rx']
            kp_soc = param['kp_soc']
            vaspy.MakeInputs.all(dst_path, incar, kp_rx, kp_soc)

    def make_list_run(self, run_file):
        """
        next.py用のlist_run.txtを作成
        """
        lines = ""
        for param in self.series:
            path = param['path']
            run_file = os.path.abspath(run_file)
            lines += "{0}    {1}\n".format(path, run_file)
        Cabinet.write_file('list_run', lines)

    def append_list_run(self, run_file):
        """
        next.py用のlist_run.txtを追加
        """
        lines = ""
        for param in self.series:
            path = param['path']
            run_file = os.path.abspath(run_file)
            lines += "{0}    {1}\n".format(path, run_file)
        Cabinet.append_file('list_run', lines)

    def mkdirs(self):
        """
        pathに記された階層ディレクトリを作成
        階層の順番はsetした順番となる
        """
        for param in self.series:
            dst_path = param['path']
            Bash.mkdir(dst_path)
