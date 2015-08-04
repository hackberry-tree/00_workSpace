#!/usr/bin/python
# -*- coding: utf-8 -*-
"""bandをplotする"""
import os
import glob
import pylab
import shutil
import unittest
import grapy
import collect_vasp
import numpy as np
#import vasp_unfold
from commopy import DataBox

PATH ="/Users/enoki/workspace/test_band/Fe3Ga_prim"

def main():
    """main"""
    polarized(PATH)

def polarized(path='.', procar='PROCAR_band', doscar='DOSCAR_polarized'):
    """polarizedのPROCARを読み込む"""
    dos = collect_vasp.Doscar(os.path.join(path, doscar))
    dos.get_data()
    e_fermi = dos.fermi_energy  # DOSCARからのEf

    band = collect_vasp.Procar(os.path.join(path, procar))
    band['energy'] = band['energy'] - e_fermi
    band.output_keys = ['kpoint_id', 'energy', "spin"]
    band.set_spin_direction()
    band.data = band.trim_data(band['spin'], 'down')

    #kp_label = band.set_turning_kpoints()
    kp_label = [[0], [0]]

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

    #ax1.scatter(band['kpoint_id'], band['energy'], s=1, c='gray',
    #            linewidths=0)
    orbitals = DataBox(band['orbitals_tot'])
    size = [float(x) * 10 for x in orbitals['dxy']]
    band_list = np.array([[band['kpoint_id'][i], band['energy'][i], size[i]]
                          for i in range(len(size)) if size[i] > 1e-2])
    ax1.scatter(band_list[:,0], band_list[:,1], s=band_list[:,2], c='g',
                 linewidths=0)

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


if __name__ == '__main__':
    main()
