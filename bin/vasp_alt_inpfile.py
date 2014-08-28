#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
何かのファイルを変更
"""
import vaspy

def main():
    """main"""
    alt_kpoints(0.1222, 'KPOINTS_soc')

def alt_kpoints(d_p, output):
    """alt KPOINTS file"""
    poscar = vaspy.Poscar('POSCAR')
    kpoints = vaspy.Kpoints(poscar.get_lattice_length(), d_p)
    kpoints.write_kpoints(output)

if __name__ == '__main__':
    main()
