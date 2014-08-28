#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Concerned with .cif file module."""
import os
import glob
import shutil
from commopy import Cabinet, Bash


class Cif2Cell(object):
    """
    Prepare POSCAR from cif files using cif2cell code.
    """

    def __init__(self, path='.'):
        """
        Find cif files at current directory and cifs/ directory.
        And initialize new_menber attribute.
        """
        self.cif_all_list = self.find_cif_files(path)
        self.new_member = []

    @staticmethod
    def conv2poscar(cif_file, out_path):
        """
        A cif file convert to POSCAR
        """
        cif_file = os.path.normpath(cif_file)
        cmd = 'cif2cell -p vasp -f {0} -o {1}/POSCAR --no-reduce \
        --vasp-print-species --vasp-cartesian-lattice-vectors \
        --vasp-format=5'.format(cif_file, out_path)
        Bash.execute(cmd)

    @staticmethod
    def find_cif_files(path):
        """
        This method find **.cif or **.cif.txt in path and path/cifs directory.
        All cif files are moved in path/cifs/ dirctory
        """
        cifs_dir = os.path.join(path, 'cifs')
        Bash.mkdir(cifs_dir)
        search_cif = glob.glob(os.path.join(path, "*.cif"))
        search_ciftxt = glob.glob(os.path.join(path, "*.cif.txt"))
        search_cwd = search_cif + search_ciftxt
        for cif_file in search_cwd:
            Bash.move(cif_file, cifs_dir)
        search_cif = glob.glob(os.path.join(cifs_dir, "*.cif"))
        search_ciftxt = glob.glob(os.path.join(cifs_dir, "*.cif.txt"))
        search_cifs = search_cif + search_ciftxt
        if len(search_cifs) == 0:
            print('I cannot find any .cif file in {0} and {1}'
                  .format(path, cifs_dir))
            exit()
        return search_cifs

    @staticmethod
    def get_name_cif(cif_file):
        """
        The name of cif file is read. This will use directory name
        in which vasp initial files are prepared.
        """
        filename = cif_file.split('/')[-1]
        fname = filename.split('.cif')[0]
        return fname

    @classmethod
    def prep_one(cls, cif_file, dst_dir):
        """
        The dst_dir is current directory,
        this method search POSCAR and reserve it.
        If not, it will search whethre dst_dir already exist or not.
        If it already exist, this module will do nothing.
        If not, new directory make and prepare POSCAR file.
        """
        if dst_dir == '.':
            Cabinet.reserve_file('POSCAR')
        elif glob.glob(dst_dir):
            line = ("\'{0}\' directory is already exist.\n"
                    "Do Nothing...\n".format(dst_dir.split('/')[-1]))
            print(line)
            return
        else:
            Bash.mkdir(dst_dir)
        cls.conv2poscar(cif_file, dst_dir)
        shutil.copy(cif_file, dst_dir)

    def prep_all(self, dst_path):
        """
        This is wrapper of prep_one.
        """
        files = self.cif_all_list
        for cif in files:
            dirname = self.get_name_cif(cif)
            self.prep_one(cif, os.path.join(dst_path, dirname))
            self.new_member.append(dirname)
