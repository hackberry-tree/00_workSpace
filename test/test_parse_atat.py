#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for parse_atat"""


import os
import glob
import shutil
import unittest

import parse_atat


__date__ = "Sep 3 2014"

TEST_PATH = ("/Users/enoki/Researches/Analysis/Codes/01_testRun")

def main():
    unittest.main()

class Test(unittest.TestCase):
    """
    テスト
    """
    PATH = os.path.join(TEST_PATH, 'atat')

    def test_get_volume(self):
        """
        体積と磁気モーメントの情報も読み取りたい
        """
        dirc = os.path.join(self.PATH, 'Fe-Ni_02')
        four_sub = [0, 1, 3, 27, 28]
        energies = os.path.join(dirc, 'energies.txt')
        inc, exc = parse_atat.get_ids_from_energies(energies, four_sub)
        print(exc)

        dirc = os.path.join(self.PATH, 'Fe-Ni_02')
        analysis = parse_atat.Analysis.from_dirc(dirc)


    def _test_strout(self):
        """
        StrOutのテスト
        """
        src = os.path.join(self.PATH, 'Fe-Ni', '0', 'str_relax.out')
        strout = parse_atat.StrOut.from_file(src)
        #print(strout.structure)
        #print(dir(strout.prim_cif))
        #strout.prim_cif
        dst = os.path.join(self.PATH, 'Fe-Ni', '0', 'str.cif')
        strout.prim_cif(dst)

    def _test_analysis(self):
        dirc = os.path.join(self.PATH, 'Fe-Ni_02')
        four_sub = [0, 1, 3, 27, 28]
        energies = os.path.join(dirc, 'energies.txt')
        inc, exc = parse_atat.get_ids_from_energies(energies, four_sub)
        print(exc)

        dirc = os.path.join(self.PATH, 'Fe-Ni_02')
        analysis = parse_atat.Analysis.from_dirc(dirc)
        path = os.path.join(dirc, 'reproduce_energy.txt')
        analysis.set_cvm_enthalpy(path)

        pre = r"""
\documentclass[12pt,a4paper]{report}
\usepackage[dvipdfmx]{graphicx}
\usepackage{longtable}
%\usepackage{amsmath}
\begin{document}
{
\begin{center}
\small
\begin{longtable}{cccc|c}
\caption{Arrangement of elements and formation enthalpy.}\\
Formula & Space Grope & lattice & sites & $\Delta H$ \\
&  & parameters & & ($\Delta H_{CVM}$)\\
\hline\hline
\endhead
"""[1:]
        end = r"""
\hline
\end{longtable}
\end{center}

\end{document}
"""[1:]
        lines = pre
        lines += analysis.to_tex_form(form_key=['Fe', 'Ni'].index, id_order=inc)
        lines += end
        with open(os.path.join(dirc, 'inc.tex'), 'w') as wfile:
            wfile.write(lines)


        lines = pre
        lines = analysis.to_tex_form(form_key=['Fe', 'Ni'].index, id_order=exc)
        lines += end
        with open(os.path.join(dirc, 'exc.tex'), 'w') as wfile:
            wfile.write(lines)



def clean_prev(path, files):
    """
    filesを消去
    """
    trush_list = glob.glob(os.path.join(path, files))
    for trush in trush_list:
        fname = os.path.join(path, trush)
        os.remove(fname)
        print("{0} is removed.".format(fname))

def clean_prev_dir(path, dirc):
    """
    dirctoryを消去
    """
    trush_list = glob.glob(os.path.join(path, dirc))
    for trush in trush_list:
        shutil.rmtree(trush)
        print("{0} is removed.".format(trush))


if __name__ == '__main__':
    main()
