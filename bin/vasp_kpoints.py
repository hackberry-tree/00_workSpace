#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
KPOINTSを取り扱う
"""
from __future__ import division

import numpy as np
from math import ceil
import argparse
from pymatgen.io.vaspio import Kpoints, Poscar
from pymatgen.symmetry.bandstructure import HighSymmKpath

from commopy import Cabinet

__date__ = "Oct 30 2014"


def main():
    """
    use argparse
    """
    parser = argparse.ArgumentParser(description='convert POSCAR')

    # convert to primitive
    #parser.add_argument('--band', nargs='?', dest='run', const=['band_path'],
    #                    action='append_const', default=[])
    parser.add_argument('--band', nargs=1, dest='run', default=[band_path],
                        action='append')

    args = parser.parse_args()

    print(args.run)

    args.run[0](int(args.run[1][0]))

def band_path(line_density=40):
    """
    for band calculation
    """
    poscar = Poscar.from_file('POSCAR', check_for_POTCAR=False)
    highsymmkp = MyHighSymmKpath(poscar.structure)
    kpts = highsymmkp.get_kpoints(line_density)
    args = {'comment': "Kpoints for band calc",
            'kpts': kpts[0],
            'num_kpts': len(kpts[0]),
            'labels': kpts[1],
            'style': 'Reciprocal',
            'kpts_weights': [1]*len(kpts[0])}
    kpoints = Kpoints(**args)
    kpoints.write_file('KPOINTS_band')

class MyHighSymmKpath(HighSymmKpath):
    def __init__(self, structure):
        HighSymmKpath.__init__(self, structure)

    def get_kpoints(self, line_density=20):
        """
        Returns:
            the kpoints along the paths in cartesian coordinates
            together with the labels for symmetry points -Wei
        """
        list_k_points = []
        sym_point_labels = []
        for b in self.kpath['path']:
            for i in range(1, len(b)):
                start = np.array(self.kpath['kpoints'][b[i - 1]])
                end = np.array(self.kpath['kpoints'][b[i]])
                # 逆格子 catesian での距離を出す
                distance = np.linalg.norm(
                    self._prim_rec.get_cartesian_coords(start) -
                    self._prim_rec.get_cartesian_coords(end))
                # 距離に比例する整数
                nb = int(ceil(distance * line_density))
                sym_point_labels.extend([b[i - 1]] + [''] * (nb - 1) + [b[i]])
                # list_k_points.extend(
                #     [self._prim_rec.get_cartesian_coords(start)
                #      + float(i) / float(nb) *
                #      (self._prim_rec.get_cartesian_coords(end)
                #       - self._prim_rec.get_cartesian_coords(start))
                #      for i in range(0, nb + 1)])
                list_k_points.extend(
                    [start + float(i) / float(nb) * (end - start)
                     for i in range(0, nb + 1)])
        return list_k_points, sym_point_labels



if __name__ == '__main__':
    main()


