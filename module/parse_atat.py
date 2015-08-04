#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ATAT関連
pymatgen の structure class を生成させる
Class:
    StrOut:
        str.out を読んで空間群などを output する
    Analysis:
        maps を実行した dirc を指定して、結果を収集する
        磁気モーメントも読めるようにする
method:
    get_ids_from_energies:
        cvm 計算との data 比較のための method
        energies.txt を読んで 構造 ID を習得
        * で comment out した構造と別けて return する

"""

import os
import glob
import numpy as np

from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__date__ = "Nov 10 2014"

def main():
    pass


class StrOut(object):
    """
    ATATのdir中にあるstr.outのparser
    pymatgenのstructureをattributeにする
    """
    def __init__(self, structure):
        self.structure = structure

    @classmethod
    def from_file(cls, src='str.out'):
        """
        ファイルからobjectを生成
        一旦poscarに変換してcartecian poscarをpymatgenで読む
        ToDo:cartecianをdirectに変換することで直接読めるようにする
        """
        posc = cls.strout2poscar(src)
        return cls(posc.structure)

    def prim_cif(self, dst):
        """
        primitive cellでのcifフォーマットをgetする
        """

        finder = SpacegroupAnalyzer(self.structure)
        structure = finder.get_primitive_standard_structure()
        structure = finder.get_conventional_standard_structure()
        cif = CifWriter(structure, find_spacegroup=True, symprec=0.1)
        cif.write_file(dst)

    @staticmethod
    def strout2poscar(src='str.out'):
        """
        str.outをposcarの記述に変換し, Poscar objをreturn
        """
        with open(src, 'r') as rfile:
            lines = rfile.readlines()
        latt_tmp = [[float(x) for x in y.split()] for y in lines[0:3]]
        trans = [[float(x) for x in y.split()] for y in lines[3:6]]
        sites_tmp = [[float(x) for x in y.split()[0:3]] for y in lines[6:]]
        elements = [x.split()[3].split('+')[0].split('-')[0] for x in lines[6:]]
        num_atoms = [elements.count(e)
                     for e in sorted(set(elements), key=elements.index)]
        latt = np.dot(np.array(trans), np.array(latt_tmp))
        sites = np.dot(np.array(sites_tmp), np.array(latt_tmp))

        posc_str = "posc_orig\n"
        posc_str += "1.00\n"
        posc_str += "\n".join(["  ".join([str(f) for f in l]) for l in latt])
        posc_str += "\n"
        posc_str += "  ".join(sorted(set(elements), key=elements.index)) + "\n"
        posc_str += "  ".join([str(d) for d in num_atoms]) + "\n"
        posc_str += "Cartesian\n"
        posc_str += "\n".join(["  ".join([str(f) for f in s]) for s in sites])

        return Poscar.from_string(posc_str)


class Analysis(object):
    """
    mapsを実行したdircを指定する
    結果をcollectして解析
    """
    def __init__(self, path_list):
        self.path_list = path_list
        self.data = self.get_data()
        self.data.sort(key=lambda x: int(x['str_id']))
        self._set_enthalpy()


    @classmethod
    def from_dirc(cls, dirc):
        """
        errorのないdirectoryのpathのlistからobjを生成
        """
        src = os.path.join(dirc, '*', 'energy')
        path_list = [os.path.dirname(p) for p in glob.glob(src)]
        err = os.path.join(dirc, '*', 'error')
        err_list = [os.path.dirname(p) for p in glob.glob(err)]
        for err in err_list:
            try:
                path_list.remove(err)
            except ValueError:
                pass
        # path_list = [os.path.join(dirc, '0')] #  for test
        return cls(path_list)


    @staticmethod
    def _get_data_from_single_dirc(dirc, src_str="str_relax.out",
                                   src_ene='energy'):
        """
        指定したdircから構造とエネルギーを読み取る
        """
        src = os.path.join(dirc, src_str)
        strout = StrOut.from_file(src)

        src = os.path.join(dirc, src_ene)
        with open(src, 'r') as rfile:
            lines = rfile.readlines()
        num_atoms = sum(strout.structure.composition.
                        to_data_dict['unit_cell_composition'].values())
        energy = float(lines[0]) / num_atoms

        analyzer = SpacegroupAnalyzer(strout.structure)
        #std_prim = analyzer.get_primitive_standard_structure()
        std_str = analyzer.get_conventional_standard_structure()
        analyzer = SpacegroupAnalyzer(std_str)
        wyckoffs = analyzer.get_symmetry_dataset()['wyckoffs']
        formula = std_str.composition.to_data_dict['unit_cell_composition']

        symbol_spg = analyzer.get_spacegroup_symbol()
        num_spg = analyzer.get_spacegroup_number()
        spg = [symbol_spg, num_spg]

        lattice = std_str.as_dict()['lattice']

        equiv_sites = analyzer.get_symmetrized_structure().equivalent_sites
        equiv_indices = analyzer.get_symmetrized_structure().equivalent_indices
        # Wycoffs labelと組み合わせたsites_groupのlistを作る
        sites_and_wyckoffs = []
        for eq_s, eq_i in zip(equiv_sites, equiv_indices):
            sites_and_wyckoffs.append({'wyckoffs': wyckoffs[eq_i[0]],
                                       'site_grp': eq_s})
        # check
            for i in range(len(eq_i)-1):
                if wyckoffs[eq_i[i]] != wyckoffs[eq_i[i+1]] or \
                           len(eq_s) != len(eq_i):
                    print("wyckoffs label is wrong !!!!")
                    print(wyckoffs)
                    print(eq_i)
                    print(len(eq_s))
                    print(dirc)
                    exit()
        return {'formula': formula, 'lattice': lattice, 'spg': spg,
                'sites_and_wyckoffs': sites_and_wyckoffs, 'energy': energy,
                'str_id': os.path.basename(dirc)}

    def get_data(self):
        """
        self.path_list中の全てのデータを収集する
        """
        print("Total valid structures are {0}".format(len(self.path_list)))
        return [self._get_data_from_single_dirc(d) for d in self.path_list]

    def _set_enthalpy(self, ref_id=[0, 1]):
        """
        self.dataにenthalpyをsetをする
        refデータが単体のみ対応
        """
        ref_data = {list(x['formula'].keys())[0]: x['energy']
                    for x in self.data if int(x['str_id']) in ref_id}
        for data in self.data:
            ref = (sum([ref_data[x]*data['formula'][x]
                       for x in data['formula'].keys()]) /
                   sum(data['formula'].values()))
            enthalpy = data['energy'] - ref
            data.update({'enthalpy': enthalpy * 1000})

    @property
    def as_dict_str_id(self):
        """
        str_idをkeyにしたdictをreturn
        """
        return {int(x['str_id']): x for x in self.data}

    def to_tex_form(self, form=0, form_key=None, id_order=None):
        """
        Returns:
            lines: "str" texのフォーマットでの構造の情報
        Args:
            form: "int" texのフォーマットを指定
            form_key: "function"
                      原子の順番をソートする時に用いる関数
                      e.g. ['Fe', 'Ni'].index
            id_order: "1d list"
                      IDの順番を指定するlist
        """
        func = [self._tex_form_00, self._tex_form_01][form]
        lines = ""
        if id_order:
            for str_id in id_order:
                data = self.as_dict_str_id[str_id]
                lines += func(data, form_key)
        else:
            for data in self.data:
                # lines += data['str_id'] + "\n" #  for test
                lines += func(data, form_key)
        return lines

    def _tex_form_00(self, data, form_key=None):
        """
        formula & spg & lattice & sites & enthalpyを出力
        """
        form_line = [self.str_formula(data['formula'], key=form_key, tex=True)]
        spg_line = [self.str_spg(data['spg'], tex=True)]
        latt_lines = self.str_lattice_list(data['lattice'], tex=True)
        sites_lines = self.str_sites_list(data['sites_and_wyckoffs'])
        try:
            enthalpy = ["{0:.1f}".format(data['enthalpy']),
                        "({0:.1f})".format(data['cvm_enthalpy'])]
        except KeyError:
            enthalpy = ["{0:.1f}".format(data['enthalpy'])]

        lenmax = max([len(latt_lines), len(sites_lines), len(enthalpy)])

        tex_lines = ""
        for i in range(lenmax):
            tex_lines += self._get_tex_table_cell(form_line, i) + " & "
            tex_lines += self._get_tex_table_cell(spg_line, i) + " & "
            tex_lines += self._get_tex_table_cell(latt_lines, i)  + " & "
            tex_lines += self._get_tex_table_cell(sites_lines, i)  + " & "
            tex_lines += self._get_tex_table_cell(enthalpy, i)  + "\\\\"
            tex_lines += "\n"
        tex_lines += "\\hline\n"
        return(tex_lines)

    def _tex_form_01(self, data, form_key=None):
        """
        formula (ID) & spg & lattice & sites & enthalpyを出力
        """
        form_line = [self.str_formula(
            data['formula'], key=form_key, tex=True),
                     "(No. " + data['str_id'] + ")"]
        spg_line = [self.str_spg(data['spg'], tex=True)]
        latt_lines = self.str_lattice_list(data['lattice'], tex=True)
        sites_lines = self.str_sites_list(data['sites_and_wyckoffs'])
        try:
            enthalpy = ["{0:.1f}".format(data['enthalpy']),
                        "({0:.1f})".format(data['cvm_enthalpy'])]
        except KeyError:
            enthalpy = ["{0:.1f}".format(data['enthalpy'])]

        lenmax = max([len(form_line), len(latt_lines), len(sites_lines),
                      len(enthalpy)])

        tex_lines = ""
        for i in range(lenmax):
            tex_lines += self._get_tex_table_cell(form_line, i) + " & "
            tex_lines += self._get_tex_table_cell(spg_line, i) + " & "
            tex_lines += self._get_tex_table_cell(latt_lines, i)  + " & "
            tex_lines += self._get_tex_table_cell(sites_lines, i)  + " & "
            tex_lines += self._get_tex_table_cell(enthalpy, i)  + "\\\\"
            tex_lines += "\n"
        tex_lines += "\\hline\n"
        return(tex_lines)


    @staticmethod
    def _get_tex_table_cell(line_list, index):
        """
        texのtableのcellを作成する為のmethod
        indexを超える場合, 空のcell("")をreturn
        """
        try:
            return line_list[index]
        except IndexError:
            return ""

    @staticmethod
    def str_formula(formula, key=None, tex=False):
        """
        formula dictをstrに変換
        keyで元素の順番を指定してソート可能
        tex=Trueの場合, tex形式で出力
        """
        if key:
            form_list = [(x, formula[x]) for x in sorted(formula, key=key)]
        else:
            form_list = formula.items()
        lines = ""
        for elem in form_list:
            if elem[1] == 1:
                lines += elem[0]
            else:
                lines += elem[0] + {True: "$_{", False: ""}[tex] + \
                         str(int(elem[1])) + {True: "}$", False: ""}[tex]
        return lines

    def set_cvm_enthalpy(self, fname='log.txt'):
        """
        CVMのlog.txtのreproduced_energyを取り込む
        ToDo:log.txtから読めるようにする
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        rep_ene = {x.split()[0].split('*')[0][3:]: x.split()[2] for x in lines}
        for data in self.data:
            data.update({'cvm_enthalpy': float(rep_ene[data['str_id']])})

    @staticmethod
    def str_spg(spg, tex=False):
        """
        spgをstr形式に変換 "spg_symbol (spg_no)"
        tex=Trueの場合spg_symbolをtex形式で出力
        """
        s = spg[0]
        i = spg[1]
        if tex:
            if '-' in s:
                pos = s.index('-')
                s = s[0:pos] + "$\\bar{"+ "{0}".format(s[pos+1]) + "}$" + \
                    s[pos+2:]
            if '_' in s:
                pos = s.index('_')
                s = s[0:pos] + "$_"+ "{0}".format(s[pos+1]) + "$" + s[pos+2:]
        return s + " ({0})".format(i)

    @staticmethod
    def str_lattice_list(lattice, tex=False):
        """
        strに変換した要素を持つlistをreturn
        tex=Trueの場合spg_symbolをtex形式で出力
        a軸と同じ場合、他の軸は省略する
        また90°の角度も省略
        """
        if lattice['a'] == lattice['b'] and lattice['a'] == lattice['c']:
            latt = {True: ["$a$=$b$=$c$={0:.3f}".format(lattice['a'])],
                    False: ["a=b=c={0:.3f}".format(lattice['a'])]}[tex]
        elif lattice['a'] == lattice['b']:
            latt = {True: ["$a$=$b$={0:.3f}".format(lattice['a']),
                           "$c$={0:.3f}".format(lattice['c'])],
                    False: ["a=b={0:.3f}".format(lattice['a']),
                            "c={0:.3f}".format(lattice['c'])]}[tex]
        else:
            latt = {True: ["$a$={0:.3f}".format(lattice['a']),
                           "$b$={0:.3f}".format(lattice['b']),
                           "$c$={0:.3f}".format(lattice['c'])],
                    False: ["$a$={0:.3f}".format(lattice['a']),
                            "$b$={0:.3f}".format(lattice['b']),
                            "$c$={0:.3f}".format(lattice['c'])]}[tex]

        labels = {True: ["$\\alpha$=", "$\\beta$=", "$\\gamma$="],
                  False: ["alpha=", "beta=", "gamma="]}[tex]
        angles = [lattice['alpha'], lattice['beta'], lattice['gamma']]
        ang = []
        for label, val in zip(labels, angles):
            if (val - 90) ** 2 > 1e-5:
                ang.append(label + "{0:.1f}".format(val))
        return latt + ang

    @staticmethod
    def str_sites_list(sites_and_wyckoffs):
        """
        sitesをstrに変換
        そのままの形式でtexに持ち込める
        """
        sites = []
        for group in sites_and_wyckoffs:
            site = group['site_grp'][0]
            for sp, _ in site.species_and_occu.items():
                line = str(sp)
                line += " (" + str(len(group['site_grp'])) + ""
                line += group['wyckoffs'] + ") "
                line += ("{0:.3f}".format(site.a)) + " "
                line += ("{0:.3f}".format(site.b)) + " "
                line += ("{0:.3f}".format(site.c))
                sites.append(line)
        return sites


def get_ids_from_energies(path='energies.txt', priority_id=[]):
    """
    構造idを並び替えたlistをreturnする
    inputはcvmのenergies.txt
    構造はidが若い順で出力
    また*でコメントアウトしている構造もreturnする
    priority_id(list)を指定しておくとそれを先頭にlistする
    attribute
        path: energies.txt file (for cvm)
        priority_id: list of IDs
    return
        included: list of IDs
        excluded: list of IDs
    ToDo: スプリットキーの'FCC'を修正
    """
    with open(path, 'r') as rfile:
        lines = rfile.readlines()
    excluded = [int(x.split()[0].split('*')[-2].split('FCC')[-1])
                 for x in lines if x[0] == '*']
    included = [int(x.split()[0].split('*')[-2].split('FCC')[-1])
                 for x in lines if x[0] != '*']
    excluded.sort()
    included.sort()
    for i in priority_id:
        included.remove(i)
    included = priority_id + included

    return included, excluded

if __name__ == "__main__":
    main()
