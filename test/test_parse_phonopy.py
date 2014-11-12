#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test for itx.py"""


import os
import glob
import shutil
import unittest
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
pylab.rcParams['animation.ffmpeg_path'] = '/opt/local/bin/ffmpeg'
from parse_phonopy import Gibbs, Composition, GibbsWithComp, Enthalpy, \
    References, MurnaghanEoS, MurnaghanWithComp

from convex_hull import FindGS
import convex_hull

__date__ = "Sep 3 2014"

TEST_PATH = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
             "99_python/01_testRun/")

def main():
    unittest.main()

class Test(unittest.TestCase):
    """
    テスト
    """
    path = os.path.join(TEST_PATH, 'phonon')

    def est_gibbs(self):
        """
        gibbsをテスト
        1. ファイルを読む
        2. atT()
        3. set_gibbs_per_atoms()
        """
        gibbs = Gibbs.from_file(os.path.join(self.path, 'Fe',
                                             'gibbs-temperature.dat'))
        gibbs.atT(1600)

        gibbs.alt_gibbs_per_atoms(10)
        print(gibbs)

    def est_comp(self):
        """
        compositionをテスト
        1. POSCARを読む
        2. fractionで出力
        3. len()
        """
        comp = Composition.from_poscar(os.path.join(self.path, 'FeNi',
                                                    'POSCAR'))
        print(comp.dict)
        print(comp['Fe'])
        print(comp['Si'])

        comp.fract = True
        print(comp['Fe'])

        print(len(comp))

    def est_murnaghanwithcomp(self):
        mwc = MurnaghanWithComp.from_directory(os.path.join(self.path, 'FeNi3'))
        print('FeNi3')
        mwc.murnaghan.alt_atP(10)
        print(mwc.murnaghan)
        print(mwc.GatT_with_comp(300, ['Fe', 'Ni', 'Si']))

    def est_ref(self):
        dir_list = [os.path.join(self.path, 'Fe')]
        dir_list.append(os.path.join(self.path, 'Ni'))
        ref = References.from_directories(dir_list)
        ref_g = ref.at_comp({'Fe':1, 'Ni':0})
        print(ref_g)


    def est_enthalpy(self):
        ref_list = [os.path.join(self.path, 'Fe')]
        #FeNi3 = GibbsWithComp.from_dirc(os.path.join(self.path, 'Fe3Ni_L12'))
        ref_list.append(os.path.join(self.path, 'Ni'))
        ref_list.append(os.path.join(self.path, 'Si'))
        compounds = glob.glob(os.path.join(self.path, '*', 'gibbs-temperature.dat'))
        compounds = [os.path.dirname(x) for x in compounds]
        print(compounds)
        refs = References.from_directories(ref_list)

        for compound in compounds:
            comp = MurnaghanWithComp.from_directory(compound)
            print(os.path.basename(compound))
            ent = Enthalpy(comp, refs)
            for p in range(0, 17, 2):
                ent['enthalpy'] = ent.enthalpy(p)
                print('{0}, {1}'.format(p, ent.data[0]['enthalpy']))

        print(ent.atT_with_comp(100, ["Fe", "Ni", "Si"]))

        #ref_FeNi3 = ref.at_comp(FeNi3.comp.dict)
        #print(ref)
        #ent = Enthalpy(FeNi3, ref)
        #print(ent)

    def test_convex(self):
        """
        アニメーション
        """
        base_list = [os.path.join(self.path, 'Fe')]
        #FeNi3 = GibbsWithComp.from_dirc(os.path.join(self.path, 'Fe3Ni_L12'))
        base_list.append(os.path.join(self.path, 'Ni'))
        base_list.append(os.path.join(self.path, 'Si'))
        ref = References.from_directories(base_list)

        compounds = glob.glob(os.path.join(self.path, '*',
                                           'gibbs-temperature.dat'))
        compounds = [os.path.dirname(x) for x in compounds]


        fig = pylab.figure()
        ax = Axes3D(fig)
        i = 0
        elems = ['Fe', 'Ni', 'Si']
        pressure = 0
        for i in [0, 2]:
            temp_list = [0, 1000, 0]
            temp = temp_list[i]
            pres_list = [0, 0, 10]
            pressure = pres_list[i]
            bases = []
            for base in base_list:
                mwc = MurnaghanWithComp.from_directory(base)
                ent = Enthalpy(mwc, ref)
                ent['enthalpy'] = ent.enthalpy(pressure)
                bases.append(ent.atT_with_comp(temp, elems))
                bases[-1][-1] *= 96.485344520851

            not_bases = []
            for compound in compounds:
                mwc = MurnaghanWithComp.from_directory(compound)
                ent = Enthalpy(mwc, ref)
                ent['enthalpy'] = ent.enthalpy(pressure)
                not_bases.append(ent.atT_with_comp(temp, elems))
                not_bases[-1][-1] *= 96.485344520851

            mwc = MurnaghanWithComp.from_directory(os.path.join(self.path,
                                                                'Fe2NiSi_5to10'))
            ent = Enthalpy(mwc, ref)
            ent['enthalpy'] = ent.enthalpy(pressure)
            heusler = ent.atT_with_comp(temp, elems)
            heusler[-1] *= 96.485344520851
            not_bases.remove(bases[0])
            not_bases.remove(bases[1])
            not_bases.remove(bases[2])

            tri = FindGS.collect_base_triangles(bases, not_bases)
            fromgs = FindGS.fromGround(tri, [heusler])
            print('{0} {1} {2}'.format(temp, fromgs[0][2], heusler[-1]))

            color = ['blue', 'magenta', 'green']
            convex_hull.draw_convex_hull(ax, bases, not_bases,
                                         ['Fe', 'Ni', 'Si'], [-80 ,5], color=color[i])
            i += 1

        def animate(i):
            ax.view_init(30, 2*i)
            fig.canvas.draw()
        FFwriter = animation.FFMpegWriter(fps=15, metadata=dict(artist='Me'),
                                          bitrate=1800)
        anim = animation.FuncAnimation(fig, animate, frames=180, interval=1,
                                       blit=False)
        anim.save('test.mp4', writer=FFwriter)

        pylab.show()

    def est_murnaghan(self):
        cev = MurnaghanEoS.from_file(os.path.join(self.path, 'Fe2NiSi_5to10',
                                             'helmholtz-volume.dat'))
        cev.alt_per_atom(8)
        cev.alt_atP(2)
        cev['G_P2'] = cev['gibbs']
        cev.alt_atP(4)

        cev['G_P4'] = cev['gibbs']
        cev.alt_atP(6)

        cev['G_P6'] = cev['gibbs']
        cev.alt_atP(8)

        cev['G_P8'] = cev['gibbs']
        cev.alt_atP(10)

        cev['G_P10'] = cev['gibbs']
        cev.output_keys = ['temp', 'G_P0', 'G_P2', 'G_P4', 'G_P6', 'G_P8', 'G_P10']
        print(cev)

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
