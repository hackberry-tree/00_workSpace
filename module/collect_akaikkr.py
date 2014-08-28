#!/opt/local/bin/python
# -*- coding: utf-8 -*-
"""
energy, c/a, mag, volumeなどの情報を読んでdata(array)を作成する
collect_vasp.pyのEnergy classにフォーマットを揃える
"""
import re
from commopy import Cabinet, DataBox


class Energy(DataBox):
    """
    path_list中のoutput fileを読んで
    energyなどのデータを収集
    collect_vaspと違ってdirc_listではなくoutputまでを含んだpath_listを使う
    """
    def __init__(self, path_list):
        self.path_list = path_list
        DataBox.__init__(self, self.get_data())
        self.output_keys = []

    def get_data(self):
        """
        self.path_list中のデータを全修得
        """
        data_box = []
        for path in self.path_list:
            result = self.get_data_single(path)
            if result:
                data_box.append(result)
        return data_box

    @staticmethod
    def get_data_single(fname):
        """
        一つのoutputファイルからデータ修得
        変数が多いので分割することも要検討
        """
        lines = Cabinet.read_file(fname)
        key_energy = r"\s*total energy=\s*([-\d\.]*)\s*"
        meta_energy = re.compile(key_energy)
        key_atoms = (r"\s*ntyp=\s*([\d]+)\s+natm=\s*([\d]+)\s+"
                     r"ncmpx=\s*([\d]+)\s*")
        meta_atoms = re.compile(key_atoms)
        key_latt = (r"\s*brvtyp=\s*([^\s]+)\s+a=\s*([\d\.]+)\s+"
                    r"c/a=\s*([\d\.]+)\s+b/a=\s*([\d\.]+)\s*")
        meta_latt = re.compile(key_latt)
        for line in lines:
            if meta_latt.match(line):
                latt_a = meta_latt.match(line).group(2)
            if meta_atoms.match(line):
                ntype = int(meta_atoms.match(line).group(1))
                num_atoms = int(meta_atoms.match(line).group(2))
            if meta_energy.match(line):
                energy = float(meta_energy.match(line).group(1))

        pos_compo = lines.index("   type of site\n")
        key_type = r"\s*type=([^\s]+)\s+.*"
        meta_type = re.compile(key_type)
        key_conc = (r"\s*component=\s*([\d]+)\s+anclr=\s*([\d\.]+)\s+"
                    r"conc=\s*([\d\.]+)\s*")
        meta_conc = re.compile(key_conc)
        i = 1
        site_id = 0
        compo = []
        while lines[pos_compo+i] != "\n":
            if meta_type.match(lines[pos_compo+i]):
                site_id += 1
                type_name = meta_type.match(lines[pos_compo+i]).group(1)
                compo.append({'site_id': site_id, 'type_name': type_name,
                              'compositions':[]})
            if meta_conc.match(lines[pos_compo+i]):
                component = meta_conc.match(lines[pos_compo+i]).group(1)
                anclr = float(meta_conc.match(lines[pos_compo+i]).group(2))
                conc = meta_conc.match(lines[pos_compo+i]).group(3)
                compo[-1]['compositions'].append({'Z': int(anclr),
                                                  'component': int(component),
                                                  'concentration': float(conc)})
            i += 1
        try:
            return {'energy': energy/num_atoms, 'num_atoms': num_atoms,
                    'latt_a': latt_a, 'compo': compo}
        except UnboundLocalError:
            print('error {0}'.format(fname))
            return None


