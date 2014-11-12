#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
KPOINTSを取り扱う
"""
from __future__ import division

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
    highsymmkp = HighSymmKpath(poscar.structure)
    kpts = highsymmkp.get_kpoints(line_density)
    args = {'comment': "Kpoints for band calc",
            'kpts': kpts[0],
            'num_kpts': len(kpts[0]),
            'labels': kpts[1],
            'style': 'Reciprocal',
            'kpts_weights': [1]*len(kpts[0])}
    kpoints = Kpoints(**args)
    kpoints.write_file('KPOINTS_band')

if __name__ == '__main__':
    main()
