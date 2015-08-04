#!/usr/bin/python
# -*- coding: utf-8 -*-
"""test"""
import os
import glob
import pylab
import shutil
import unittest
import grapy
import collect_vasp
#import vasp_unfold
from commopy import DataBox


def main():
    """main"""
    unittest.main()

TEST_PATH = '/Users/enoki/Researches/Analysis/Codes/01_testRun/'

class Dos(unittest.TestCase): #pylint: disable=R0904
    """DOSCAR用"""
    path = os.path.join(TEST_PATH, 'band')

    def dos(self):
        """dosを読む"""
        path = os.path.join(self.path, 'data')
        dos = collect_vasp.Doscar(os.path.join(path, 'DOSCAR_polarized'))
        dos.get_data()
        print(dos)


class Procar(unittest.TestCase): #pylint: disable=R0904
    """PROCAR用"""
    path = os.path.join(TEST_PATH, 'band')

    def _test_vs_vasp_unfold(self):
        path = os.path.join(self.path, 'data')
        band = collect_vasp.Procar(os.path.join(path, 'PROCAR_band'))
        print(band['kpoint'])

        band2 = vasp_unfold.parse_procar(os.path.join(path, 'PROCAR_band'))
        print(band2[5][0])



    def _test_band_polarized(self):
        """polarizedのPROCAR用"""
        path = os.path.join(self.path, 'data')

        dos = collect_vasp.Doscar(os.path.join(path, 'DOSCAR_polarized'))
        dos.get_data()
        e_fermi = dos.fermi_energy  # DOSCARからのEf

        band = collect_vasp.Procar(os.path.join(path, 'PROCAR_band'))
        band['energy'] = band['energy'] - e_fermi
        band.output_keys = ['kpoint_id', 'energy', "spin"]
        band.set_spin_direction()
        band.data = band.trim_data(band['spin'], 'down')

        kp_label = band.set_turning_kpoints()

        # plot
        plt = pylab.figure(figsize=(7, 7/1.618))
        default_font = {'font.family': 'Times New Roman'}
        pylab.rcParams.update(default_font)
        ax1 = plt.add_subplot(111)
        ax1.set_ylabel("Energy from $E_F$(meV)")
        ax1.set_xlabel("BZ direction")

        # 枠
        for i in kp_label[0]:
            pylab.axvline(x=i, ls=':', color='gray')
        pylab.axhline(y=0, ls=':', color='gray')

        ax1.scatter(band['kpoint_id'], band['energy'], s=1, c='gray',
                    linewidths=0)
        orbitals = DataBox(band['orbitals_tot'])
        size = [float(x) * 10 for x in orbitals['dxy']]
        ax1.scatter(band['kpoint_id'], band['energy'], s=size, c='g',
                    linewidths=0)
        size = [float(x) * 10 for x in orbitals['dx2']]
        ax1.scatter(band['kpoint_id'], band['energy'], s=size, c='b',
                    linewidths=0)
        ax1.set_ylim(-5, 5)
        ax1.set_xlim(0, band['kpoint_id'][-1])
        pylab.xticks(kp_label[0], kp_label[1])
        pylab.show()
        #pylab.savefig(os.path.join(path, 'band.eps'))

    def test_band_soc02(self):
        path = os.path.join(self.path, 'FeCo_new')
        path = os.path.join(self.path, 'Fe16N2/elem_Fe12Co4N2_prim')
        path = os.path.join(self.path, 'Fe16N2/Fe16N2')
        dos = collect_vasp.Doscar(os.path.join(path, 'dos'))
        dos.get_data()
        ef = dos.fermi_energy
        print(ef)

        band_001 = collect_vasp.ProcarNonC(os.path.join(path, 'band001'))
        kp_label = band_001.set_turning_kpoints_pmg(
            os.path.join(path, 'KPOINTS'))
        band_001.set_spin_direction()
        band_001.data = band_001.trim_data(band_001['spin'], 'down')
        band_001['energy'] = band_001['energy'] - ef

        # plot set style
        plt = pylab.figure(figsize=(7, 7/1.618))
        default_font = {'font.family': 'Times New Roman'}
        pylab.rcParams.update(default_font)
        ax1 = plt.add_subplot(111)
        ax1.set_ylabel("Energy from $E_F$(meV)")
        ax1.set_xlabel("BZ direction")

        # 枠
        for i in kp_label[0]:
            pylab.axvline(x=i, ls=':', color='gray')
        pylab.axhline(y=0, ls=':', color='gray')

        orbitals = DataBox(band_001['orbitals_tot'])
        mark_size = [float(x) * 10 for x in orbitals['dxy']]
        ax1.scatter(band_001['kpoint_id'], band_001['energy'], s=mark_size,
                    c='g', linewidths=0)

        mark_size = [float(x) * 10 for x in orbitals['dx2']]
        ax1.scatter(band_001['kpoint_id'], band_001['energy'], s=mark_size,
                    c='b', linewidths=0)

        band_100 = collect_vasp.ProcarNonC(os.path.join(path, 'band_nc'))
        band_100.set_spin_direction()
        band_100.data = band_100.trim_data(band_100['spin'], 'down')
        band_100['energy'] = band_100['energy'] - ef
        orbitals = DataBox(band_100['orbitals_tot'])
        mark_size = [float(x) * 10 for x in orbitals['dxy']]
        ax1.scatter(band_100['kpoint_id'], band_100['energy'], s=mark_size,
                    c='r', linewidths=0)

        mark_size = [float(x) * 10 for x in orbitals['dx2']]
        ax1.scatter(band_100['kpoint_id'], band_100['energy'], s=mark_size,
                    c='y', linewidths=0)

        ax1.set_ylim(-1, 1)
        pylab.xticks(kp_label[0], kp_label[1])
        pylab.show()

    def _test_band_sum02(self):
        """ef以下のbandのenergyのsumを取る"""
        path = os.path.join(self.path, 'FeCo_new')
        dos_001 = collect_vasp.Doscar(os.path.join(path, 'dos'))
        dos_001.get_data()
        ef_001 = dos_001.fermi_energy

        dos_100 = collect_vasp.Doscar(os.path.join(path, 'dos'))
        dos_100.get_data()
        ef_100 = dos_100.fermi_energy

        band_001 = collect_vasp.ProcarNonC(os.path.join(path, 'band001'))
        band_001.set_spin_direction()
        band_001.data = band_001.trim_data(band_001['spin'], 'down')
        band_001['energy'] = band_001['energy'] - ef_001
        sum_e001 = band_001.get_sum_energies()

        band_100 = collect_vasp.ProcarNonC(os.path.join(path, 'band100'))
        band_100.set_spin_direction()
        band_100.data = band_100.trim_data(band_100['spin'], 'down')
        band_100['energy'] = band_100['energy'] - ef_100
        sum_e100 = band_100.get_sum_energies()

        sum_e001['dif'] = sum_e100['Sum_E'] - sum_e001['Sum_E']

        plt = grapy.Vertical(1)
        #orbitals = DataBox(band_001['orbitals_tot'])
        #mark_size = [float(x) * 10 for x in orbitals['dxy']]
        plt.ax1.scatter(sum_e001['kpoint_id'], sum_e001['dif'], s=5,
                        c='g', linewidths=0)
        plt.ax1.set_ylim(-0.5, 0.5)
        plt.plot('show')


    def _test_band_soc(self):
        """soc用のband"""
        path = os.path.join(self.path, 'data_soc')
        dos_001 = collect_vasp.Doscar(os.path.join(path, 'dos_001'))
        dos_001.get_data()
        ef_001 = dos_001.fermi_energy
        print(ef_001)

        dos_100 = collect_vasp.Doscar(os.path.join(path, 'dos_100'))
        dos_100.get_data()
        ef_100 = dos_100.fermi_energy
        print(ef_100)

        band_001 = collect_vasp.ProcarNonC(os.path.join(path, 'band_001'))
        kp_label = band_001.set_turning_kpoints()
        band_001.set_spin_direction()
        band_001.data = band_001.trim_data(band_001['spin'], 'down')
        band_001['energy'] = band_001['energy'] - ef_001

        # plot set style
        plt = pylab.figure(figsize=(7, 7/1.618))
        default_font = {'font.family': 'Times New Roman'}
        pylab.rcParams.update(default_font)
        ax1 = plt.add_subplot(111)
        ax1.set_ylabel("Energy from $E_F$(meV)")
        ax1.set_xlabel("BZ direction")

        # 枠
        for i in kp_label[0]:
            pylab.axvline(x=i, ls=':', color='gray')
        pylab.axhline(y=0.03, ls=':', color='gray')

        orbitals = DataBox(band_001['orbitals_tot'])
        mark_size = [float(x) * 10 for x in orbitals['dxy']]
        ax1.scatter(band_001['kpoint_id'], band_001['energy'], s=mark_size,
                    c='g', linewidths=0)

        mark_size = [float(x) * 10 for x in orbitals['dx2']]
        ax1.scatter(band_001['kpoint_id'], band_001['energy'], s=mark_size,
                    c='b', linewidths=0)

        band_100 = collect_vasp.ProcarNonC(os.path.join(path, 'band_100'))
        band_100.set_spin_direction()
        band_100.data = band_100.trim_data(band_100['spin'], 'down')
        band_100['energy'] = band_100['energy'] - ef_100
        orbitals = DataBox(band_100['orbitals_tot'])
        mark_size = [float(x) * 10 for x in orbitals['dxy']]
        ax1.scatter(band_100['kpoint_id'], band_100['energy'], s=mark_size,
                    c='r', linewidths=0)

        mark_size = [float(x) * 10 for x in orbitals['dx2']]
        ax1.scatter(band_100['kpoint_id'], band_100['energy'], s=mark_size,
                    c='y', linewidths=0)

        ax1.set_ylim(-.5, .5)
        pylab.xticks(kp_label[0], kp_label[1])
        pylab.show()

    def _test_band_sum(self):
        """ef以下のbandのenergyのsumを取る"""
        path = os.path.join(self.path, 'data_soc')
        dos_001 = collect_vasp.Doscar(os.path.join(path, 'dos_001'))
        dos_001.get_data()
        ef_001 = dos_001.fermi_energy

        dos_100 = collect_vasp.Doscar(os.path.join(path, 'dos_100'))
        dos_100.get_data()
        ef_100 = dos_100.fermi_energy

        band_001 = collect_vasp.ProcarNonC(os.path.join(path, 'band_001'))
        band_001.set_spin_direction()
        band_001.data = band_001.trim_data(band_001['spin'], 'down')
        band_001['energy'] = band_001['energy'] - ef_001 -0.05
        sum_e001 = band_001.get_sum_energies()

        band_100 = collect_vasp.ProcarNonC(os.path.join(path, 'band_100'))
        band_100.set_spin_direction()
        band_100.data = band_100.trim_data(band_100['spin'], 'down')
        band_100['energy'] = band_100['energy'] - ef_100 -0.05
        sum_e100 = band_100.get_sum_energies()

        sum_e001['dif'] = sum_e100['Sum_E'] - sum_e001['Sum_E']

        plt = grapy.Vertical(1)
        #orbitals = DataBox(band_001['orbitals_tot'])
        #mark_size = [float(x) * 10 for x in orbitals['dxy']]
        plt.ax1.scatter(sum_e001['kpoint_id'], sum_e001['dif'], s=5,
                        c='g', linewidths=0)
        plt.ax1.set_ylim(-0.5, 0.5)
        plt.plot('show')



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
