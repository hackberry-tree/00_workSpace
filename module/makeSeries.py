#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
要整理
"""
import os
import re
import copy
import glob
import solid
from commopy import TrialRun

Combi_list = [["Co"], ["Al", "Ga", "Ge", "Si", "Sn"]]

def main():
    dirc = '.'
    Trial = TrialRun(testRun)
    Trial.whetherExist(os.path.join(dirc, 'original'))

    #Heusler()
    B2()


def testRun():
    dirc = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
            "99_python/01_testRun/makeSeriesAKAI/Heusler")
    comp_list = [["Co"], ["Al", "Ga", "Ge"]]
    conc_list = [0, 50, 100]
    latt_list = [9.0, 9.2, 9.4]
    cova_list = [0.9, 1.0, 1.1]
    Heusler(dirc, comp_list, conc_list, latt_list, cova_list)


def B2():
    os.system('mkdir Calc')
    os.system('cp -r originals Calc')
    os.chdir('Calc')
    variables_list_latt = [4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1,
                           5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8]
    Latt = AKAI_Latt('latt_dep', variables_list_latt)
    for var in Latt.variables_list:
        Latt.mkforms(var)
        output_dir = Latt.dir
        Latt.cpOriginals(Latt.original_dir, output_dir)
        Latt.convert(output_dir, var)
        Latt.addList(output_dir, 'run.sh')


def Heusler(dirc, combi_list, conc_list, latt_list, cova_list):
    #os.system('cp -r %s/originals %s/series'%(dirc, dirc))
    variables_list_latt = [9.0, 9.2, 9.4, 9.6, 9.8, 10.0,
                           10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 11.6]

    Patterns = MakePattern()
    Combi_1_2 = Patterns.combination_1_2(combi_list)
    Conc = AKAI_Conc('Conc_dep', conc_list)
    Latt = AKAI_Latt('latt_dep', latt_list)
    CovA = AKAI_CovA('CovA_dep', cova_list)
    CompY01 = AKAI_CompY01('Comp_dep')
    CompY02 = AKAI_CompY02('Comp_dep')
    CompZ01 = AKAI_CompZ01('Comp_dep')

    for list in Combi_1_2:
        combi_dir = list['pair']
        Y = list['comp1']['ele1']
        Z1 = list['comp1']['ele2']
        Z2 = list['comp2']['ele2']

        for var0 in Conc.variables_list:
            Conc.mkforms(var0)
            for var1 in Latt.variables_list:
                Latt.mkforms(var1)
                for var2 in CovA.variables_list:
                    CovA.mkforms(var2)
                    output_dir = os.path.join(dirc, 'series', combi_dir,
                                              Conc.dir, Latt.dir, CovA.dir)
                    print(output_dir)
                    Latt.cpOriginals(os.path.join(dirc, 'originals'),
                                     output_dir)
                    Conc.convert(output_dir, var0)
                    Latt.convert(output_dir, var1)
                    CovA.convert(output_dir, var2)
                    Latt.addList(output_dir, 'run.sh')
                    CompY01.convert(output_dir, Z1)
                    CompY02.convert(output_dir, Z2)
                    CompZ01.convert(output_dir, Y)


class SeriesBox:
    """
    パラメータを収納する箱
    convfunc_listには変えるべきファイルと変数とフォームを記載した組み合わせを入力
    """
    def __init__(self, title, convfunc_list, variable_list,
                 dirform_func, original_dir):
        self.title = title
        self.convfunc_list = convfunc_list
        self.variables_list = variable_list
        self.dirform_func = dirform_func
        self.original_dir = self.original_dir

    def mkforms(self, var):
        """
        postとdirnameをfuncに基づいて変換
        self.convfunc_list => [{'filename':  ,'pre':  ,
                                'post':  , 'post_func':}]
        self.dir = dir
        """
        for dict in self.convfunc_list:
            dict.update({'post': dict['post_func'](var)})
        self.dir = self.dirform_func(var)

    def mkdir(self, path):
        os.system('mkdir -p ' + path)

    def cpOriginals(self, original_dir, output_dir):
        self.mkdir(output_dir)
        os.system('cp -r %s/* %s' % (original_dir, output_dir))

    def sed(self, orig_lines, pre, post):
        new_lines = []
        for line in orig_lines:
            new_lines.append(re.sub(pre, post, line))
        return new_lines

    def convert(self, path, var):
        self.mkforms(var)
        for conv_dict in self.convfunc_list:
            lines = readFile('%s/%s' % (path, conv_dict['filename']))
            lines = self.sed(lines, conv_dict['pre'], conv_dict['post'])
            writeFile('%s/%s' % (path, conv_dict['filename']), lines)

    def addList(self, path, run_file):
        line = "%s %s" % (path, run_file)
        os.system("echo %s >> list_run.txt" % (line))


class AKAI_Latt(SeriesBox):
    #===for AKAI-KKR latt-dep=================
    def __init__(self, title, variables_list):
        self.input_file = "input"
        self.original_dir = "originals"
        convfunc_list = [{'filename': self.input_file,
                          'pre': "_Latt_", 'post_func': self.convfuncLatt}]
        #variables_list = [9.0]
        SeriesBox.__init__(self, title, convfunc_list,
                           variables_list, self.convfuncDirName,
                           self.original_dir)

    def convfuncLatt(self, var):
        post = "%s" % var
        return post

    def convfuncDirName(self, var):
        post = "latt_%.2f" % var
        return post


class AKAI_CovA(AKAI_Latt):
    #===for AKAI-KKR c/a-dep================
    def __init__(self, title, variables_list):
        self.input_file = "input"
        self.original_dir = "originals"
        convfunc_list = [{'filename': self.input_file, 'pre': "_CovA_",
                          'post_func': self.convfuncLatt}]
        SeriesBox.__init__(self, title, convfunc_list, variables_list,
                           self.convfuncDirName, self.original_dir)

    def convfuncDirName(self, var):
        post = "CovA_%.2f" % var
        return post


class AKAI_Conc(SeriesBox):
    #===for AKAI-KKR concentration-dep======
    def __init__(self, title, variables_list):
        self.input_file = "input"
        self.original_dir = "originals"
        convfunc_list = [{'filename': self.input_file, 'pre': "_Conc01_",
                          'post_func': self.convfuncConc01}]
        convfunc_list.append({'filename': self.input_file, 'pre': "_Conc02_",
                              'post_func': self.convfuncConc02})
        SeriesBox.__init__(self, title, convfunc_list, variables_list,
                           self.convfuncDirName, self.original_dir)

    def convfuncConc01(self, var):
        post = "%d" % var
        return post

    def convfuncConc02(self, var):
        var = 100 - var
        post = "%d" % var
        return post

    def convfuncDirName(self, var):
        post = "Conc_%d" % var
        return post
#========================================


class AKAI_CompY01(SeriesBox):
    #===for AKAI-KKR composition-dep=========
    def __init__(self, title):
        self.input_file = "input"
        self.original_dir = "originals"
        convfunc_list = [{'filename': self.input_file, 'pre': "_ElemY01_",
                          'post_func': self.convfuncElem}]
        convfunc_list.append({'filename': self.input_file, 'pre': "_ZY01_",
                              'post_func': self.convfuncZ})
        variables_list = []
        SeriesBox.__init__(self, title, convfunc_list, variables_list,
                           self.convfuncDirName, self.original_dir)

    def convfuncElem(self, var):
        post = "%s" % var
        return post

    def convfuncZ(self, var):
        var = Elements[var]['Z']
        post = "%d" % var
        return post

    def convfuncDirName(self, var):
        post = "%s" % var
        return post
#========================================


class AKAI_CompY02(AKAI_CompY01):
    #===for AKAI-KKR composition-dep=========
    def __init__(self, title):
        AKAI_CompY01.__init__(self, title)
        convfunc_list = [{'filename': self.input_file, 'pre': "_ElemY02_",
                          'post_func': self.convfuncElem}]
        convfunc_list.append({'filename': self.input_file, 'pre': "_ZY02_",
                              'post_func': self.convfuncZ})
        variables_list = []
        SeriesBox.__init__(self, title, convfunc_list, variables_list,
                           self.convfuncDirName, self.original_dir)
#========================================


class AKAI_CompZ01(AKAI_CompY01):
    #===for AKAI-KKR composition-dep=========
    def __init__(self, title):
        AKAI_CompY01.__init__(self, title)
        convfunc_list = [{'filename': self.input_file, 'pre': "_ElemZ01_",
                          'post_func': self.convfuncElem}]
        convfunc_list.append({'filename': self.input_file, 'pre': "_ZZ01_",
                              'post_func': self.convfuncZ})
        variables_list = []
        SeriesBox.__init__(self, title, convfunc_list, variables_list,
                           self.convfuncDirName, self.original_dir)
#========================================


class Combinatorial(object):
    def __init__(self, *combi_list):
        """
        combi_listの全ての組み合わせのリストを作成
        各要素はdict形式で'elements', 'compositon'をkeyに持つ
        """
        tmp_list = MakePattern.make_tree(*combi_list)
        self.compositions = [{'compositon': "".join(x),
                              'elements': x} for x in tmp_list]

    def set_formula(self, *num_atoms):
        """
        elementの個数をself.compostionsに設定して
        formulaを作成
        """
        for comp in self.compositions:
            comp.update({'num_atoms': list(num_atoms)})
            formula = self.make_formula(comp['elements'], num_atoms)
            comp.update({'formula': formula})

    @staticmethod
    def make_formula(elements, num_atoms):
        """
        Make chemical formula.
        """
        if len(elements) != len(num_atoms):
            print("ERROR: num_atoms list dose not match elements list !!!")
            return

        formula = ""
        for element, num in zip(elements, num_atoms):
            if num == 1:
                num = ''
            formula += "{0}{1}".format(element, num)
        return formula


class MakePattern(object):  # self.compos{'composition':(X,Y)}
    def mkDirList(self, path, list):
        #[os.path.join(path, x) for x in list]で置き換える！！
        self.dir_list = []
        for item in list:
            self.dir_list.append('%s/%s' % (path, item))

    @staticmethod
    def make_tree(*in_list):
        """
        個々の入力リストから、それぞれ一つずつ要素を取り出して
        すべての組み合わせに対するリスト(ツリー)を出力する
        例:
        a = ['a', 'b']
        b = ['c', 'd']
        c = ['e']
        print(MakePattern.make_tree(a,b,c))
        >>[['a', 'c', 'e'], ['a', 'd', 'e'], ['b', 'c', 'e'], ['b', 'd', 'e']]
        """
        dst_list = [[x] for x in in_list[0]]
        for src_list in in_list[1:]:
            dst_list = [x + [y] for x in dst_list for y in src_list]
        return dst_list

    @staticmethod
    def nCr(n, r):
        """
        1-nまでの整数から、r個の要素を取り出す場合の全ての組み合わせを出力する
        """
        out_list = []
        for i in range(0, n - r + 1):
            out_list.append([i])
        for i in range(0, r - 1):
            tmp_list = []
            for output in out_list:
                for j in range(max(output) + 1, n - r + 2 + i):
                    tmp_list.append(copy.deepcopy(output) + [j])
            out_list = copy.deepcopy(tmp_list)
        return out_list

    @classmethod
    def nCrList(cls, in_list, r):
        """
        listの中からr個選ぶ全ての組み合わせをreturn
        """
        n = len(in_list)
        index_list = cls.nCr(n, r)
        out_list = []
        for index in index_list:
            tmp_list = []
            for i in range(0, r):
                tmp_list.append(in_list[index[i]])
            out_list.append(tmp_list)
        return out_list

    def combination_1_2(self, combi_list):
        """
        [Y], [Z] => [YZ1, YZ2]
        self.pares_list[{'pair':'YZ1-YZ2',
        'comp1':{comp:'YZ1',ele1:'Y',ele2:'Z1'},
        'comp2':{comp:'YZ2',ele1:'Y',ele2:'Z2'}]
        """
        Y_list = combi_list[0]
        Z_list = self.nCrList(combi_list[1], 2)
        patterns = self.makeTree(Y_list, Z_list)

        pairs_list = []

        for pattern in patterns:
            dict = {}
            [Y, [Z1, Z2]] = pattern
            dict.update({'pair': '%s%s-%s%s' % (Y, Z1, Y, Z2)})
            dict.update({'comp1': {'comp': Y + Z1, 'ele1': Y, 'ele2': Z1}})
            dict.update({'comp2': {'comp': Y + Z2, 'ele1': Y, 'ele2': Z2}})
            pairs_list.append(dict)
        return pairs_list


def readFile(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close
    return lines


def writeFile(fname, lout):
    out = open(fname, 'w')
    for line in lout:
        out.write(line)
    out.close

#===Elements Data Base====================
Elements = solid.ELEMENTS
#==========================================

if __name__ == "__main__":
    main()
