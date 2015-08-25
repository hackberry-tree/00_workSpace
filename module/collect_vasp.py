#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
energy, c/a, mag, volumeなどの情報を読んでdata(array)を作成する
"""
import os
import re
import glob
import pylab
import numpy as np
from collections import Counter
import vaspy
import solid
from commopy import Cabinet, Array, DataBox
from fitting_analysis import FitData
from pymatgen.io.vaspio import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Mg_fit_results(DataBox):
    """
    Murnaghan fit の結果を読む
    """
    def __init__(self, data):
        DataBox.__init__(self, data)
        self.output_keys = []

    @classmethod
    def from_dir_list(cls, dir_list):
        """
        dir_list から結果を収集
        """
        data_box = []
        for dirc in dir_list:
            results = {}
            results.update(
                cls.parse_fit_results(os.path.join(dirc, 'fit_results.dat')))
            results.update(cls.parse_poscar(os.path.join(dirc, 'POSCAR')))
            results.update({'path': dirc})
            data_box.append(results)
        return cls(data_box)

    @staticmethod
    def parse_fit_results(path='fit_results.dat'):
        """
        fit_results.dat を読み込む
        """
        with open(path, 'r') as rfile:
            lines = rfile.readlines()
        energy = float(lines[0].split()[1])
        B = float(lines[1].split()[1])
        dB = float(lines[1].split()[3])
        B1 = float(lines[2].split()[1])
        dB1 = float(lines[2].split()[3])
        volume = float(lines[3].split()[1])
        return {'energy0': energy, 'energy': energy, 'B': B, 'B1': B1,
                'volume0': volume, 'volume': volume, 'dB': dB, 'dB1': dB1}

    @staticmethod
    def parse_poscar(path='POSCAR'):
        """
        POSCAR を読み込む
        """
        with open(path, 'r') as rfile:
            lines = rfile.readlines()
        elements = lines[5].split()
        num_atoms = [int(x) for x in lines[6].split()]
        return {'elements': elements, 'num_atoms': num_atoms}

    def alt_energy_at_fixed_volume(self, new_volume):
        """
        fitting 結果に基づいて 一定体積の energy に変更
        'C' 原子は数えず、ユニットセルを一定とする
        """
        for data in self.data:
            alt_volume = new_volume * (1 - data['comp_dict_f']['C'])
            new_energy = FitData.Murnaghan_func(
                [data['energy0'], data['B'], data['B1'], data['volume0']],
                alt_volume)
            data['energy'] = new_energy
            data['volume'] = alt_volume

    def set_comp_dict_f(self):
        """
        組成比を分率で作成
        """
        comp_all = set([x for y in self['elements'] for x in y])
        for data in self.data:
            comp_dict = {elem: 0 for elem in comp_all}
            sum_atoms = float(sum(data['num_atoms']))
            comp_dict_r = {elem: num/sum_atoms for elem, num
                           in zip(data['elements'], data['num_atoms'])}
            comp_dict.update(comp_dict_r)
            data.update({'comp_dict_f': comp_dict})

    def make_energies_txt(self):
        """
        構造の名前は str.out と atoms を使って作成する
        """
        with open('atoms', 'r') as rfile:
            lines = rfile.readlines()
        label = 'ABCDEFG'
        atoms = {label[i]:lines[i][:-1] for i in range(len(lines))}
        self.data.sort(key=lambda x: x['path'])

        lines_out = ""
        for data in self.data:
            with open(os.path.join(data['path'], "str.out"), 'r') as rfile:
                lines = rfile.readlines()
            elements = Counter([x.split()[-1] for x in lines[6:]])
            lines_out += data['path'] + "*"
            for key in sorted(atoms.keys()):
                lines_out += key + str(elements[atoms[key]])
            lines_out += " 1 " + str(data['energy']) + " 1\n"
        with open("energies.txt", 'w') as wfile:
            wfile.write(lines_out)

    def __str__(self):
        """
        Print self.data
        特定の要素のみを書き出したい場合は
        self.output_keysに指定する
        defaultでは全てを書き出す
        """
        keys = self.output_keys
        if not keys:
            keys = self.data[0].keys()
        types_dict = [x for x in keys if type(self.data[0][x]) is dict]
        out_lines = "\t".join(keys) + "\n"
        for key in types_dict:
            head = "".join(key.split('_dict')) + "_"
            label = [head + str(x) for x in sorted(self.data[0][key])]
            alt_label = "\t".join(label)
            out_lines = out_lines.replace(key, alt_label)
        for data in self.data:
            line = []
            for key in keys:
                if type(data[key]) is dict:
                    data[key] = "\t".join([str(data[key][x]) for x
                                           in sorted(data[key])])
                line.append(str(data[key]))
            out_lines += "\t".join(line) + "\n"
        return out_lines


class Procar(DataBox):
    """
    PROCARを読む
    """
    orbital_lists = {'s': [1, None, None], 'py': [2, None, None],
                     'pz': [3, None, None], 'px': [4, None, None],
                     'dxy': [5, 'green', r"dxy"],
                     'dyz': [6, 'salmon', r"dyz"],
                     'dz2': [7, 'purple', r"d3z^{/=7 2}-r^{/=7 2}"],
                     'dxz': [8, 'red', r"dxz"],
                     'dx2': [9, 'blue', r"dx^{/=7 2}-y^{/=7 2}"]} #  色とラベルを指定

    def __init__(self, procar, poscar=None):
        """
        self.dataの構成
        'kpoint_id', 'orbitals_phase', 'kpoint', 'energy', 'spin'
        'orbitals_tot', 'occupancy'
        'dk', 'xaxis'を追加
        'dk'は前のk点からのk空間上の距離,
        'xaxis'はそれを積算したものでband描写の横軸に用いる
        poscarを読み込ませることで逆格子cellの各軸の長さを反映させて,
        bandの横軸の長さを対応させることができる
        defaultではposcarは読まずに、a*, b*, c*の長さは全て1
        """
        with open(procar, 'r') as flow_procar:
            line = flow_procar.readline()
            line = flow_procar.readline()
            (self.num_kpoints,
             self.num_bands,
             self.num_atoms) = self.get_meta(line)
            print(self.get_meta(line))
            DataBox.__init__(self, self.get_data(flow_procar))

    def set_turning_kpoints(self):
        """
        k点のsamplingの方向が切り替わる位置を読む
        """
        kpoints = [self.data[i]['kpoint']
                   for i in range(0, len(self.data), self.num_bands)]
        kpoints = np.array(kpoints)
        delta_k = np.sign([kpoints[i+1]-kpoints[i]
                           for i in range(0, len(kpoints)-1)])
        turn_pt = [0] + [i+1 for i in range(0, len(delta_k)-1)
                         if not (delta_k[i] == delta_k[i+1]).all()]
        turn_pt += [len(kpoints) - 1]
        label = [self.get_kpoint_label(kpoints[i]) for i in turn_pt]
        return turn_pt, label

    @staticmethod
    def set_turning_kpoints_pmg(fname='KPOINTS'):
        """
        k点のsamplingの方向が切り替わる位置を読む
        pmgから生成したkpointsの場合
        KPOINTSにラベルがふられるのそこから直接読む
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        turn_pt = []
        label = []

        for i in range(len(lines)):
            if len(lines[i].split()) == 5:
                turn_pt.append(i-3)
                label.append(lines[i].split()[4])
        return turn_pt, label

    @staticmethod
    def get_kpoint_label(kpoint, structure='bct'):
        """kpointのラベルを割当てる"""
        label_bct = {0.0: {0.0: {0.0: r"$\Gamma$",
                                 0.5: "Z"},
                           0.5: {0.0: "X",
                                 0.5: "R"}},
                     0.5: {0.0: {0.0: "X",
                                 0.5: "R"},
                           0.5: {0.0: "M",
                                 0.5: "A"}}}
        label = {'bct': label_bct}
        return label[structure][kpoint[0]][kpoint[1]][kpoint[2]]

    @staticmethod
    def get_meta(line):
        """
        parameterを読み取る
        """
        meta = re.compile(r"#[^:]+:([^#]+)#[^:]+:([^#]+)#[^:]+:(.+)$")
        match_line = meta.match(line)
        num_kpoints = int(match_line.group(1))
        num_bands = int(match_line.group(2))
        num_atoms = int(match_line.group(3))
        return num_kpoints, num_bands, num_atoms

    @staticmethod
    def get_kpoint(flow_lines):
        """
        k点(h,k,l)の読む
        flow_linesはk点記載の箇所まで更新
        """
        keywords = (r"\s*k-point\s*([\d]+)*\s*:"
                    r"\s*([^\s]+)\s*([^\s]+)\s*([^\s]+).*")
        meta_kpoint_header = re.compile(keywords)
        line = flow_lines.readline()
        while not meta_kpoint_header.match(line):
            line = flow_lines.readline()  # headerまでpass
        match_line = meta_kpoint_header.match(line)
        kpoint_id = int(match_line.group(1))
        kpoint = [float(match_line.group(x)) for x in range(2, 5)]
        return kpoint_id, kpoint

    def get_1kpoint(self, flow_lines, num_bands, num_atoms):
        """
        一つのk点に対してband(energy)のセットを読み込む
        readline()関数を使うことで一行ごと処理を行う
        """
        kpoint = self.get_kpoint(flow_lines)
        meta_band_energy = re.compile(r".+energy\s+([^#]+)# occ\.\s(.+)$")
        points = []
        for _ in range(0, num_bands):  # バンドの数だけ繰り返し処理
            line = flow_lines.readline()
            while not meta_band_energy.match(line):
                line = flow_lines.readline()
            energy = float(meta_band_energy.match(line).group(1))
            occ = float(meta_band_energy.match(line).group(2))
            point = {"kpoint_id": kpoint[0], "kpoint": kpoint[1],
                     "energy": energy, "occupancy": occ}
            point.update(self.get_orbital_part(flow_lines, num_atoms))
            points.append(point)
        return points

    def get_orbital_part(self, flow_lines, num_atoms):
        """
        軌道のweightと位相を読む
        """
        orbital = self.get_weight_part(flow_lines, num_atoms)
        phase = self.get_phase_part(flow_lines, num_atoms)
        orbital.update(phase)
        return orbital

    def get_weight_part(self, flow_lines, num_atoms):
        """
        軌道のweightの部分を読む
        """
        meta_orbital_header = re.compile(r"ion\s+s.+tot$")
        line = flow_lines.readline()
        while not meta_orbital_header.match(line):
            line = flow_lines.readline()  # orbital headerまでpass
        orbitals_labels = line.split()
        orbitals = self.get_single_weight(flow_lines, num_atoms,
                                          orbitals_labels)
        return {"orbitals_ions": orbitals[0], "orbitals_tot": orbitals[1]}

    @staticmethod
    def get_single_weight(flow_lines, num_atoms, orbitals_labels):
        """
        一つのbandに対して各ionのwightとtotまでの値を読む
        ionそれぞれの値 orbitals_ions (list表記)
        すべてのionについての和 orbitals_tot
        get_weight_partで繰り返し使
        """
        orbitals_ions = {x: [] for x in orbitals_labels}
        for _ in range(0, num_atoms):
            values = flow_lines.readline().split()
            for i in range(0, len(orbitals_labels)):
                orbitals_ions[orbitals_labels[i]].append(values[i])
        line = flow_lines.readline()
        values = line.split()
        orbitals_tot = {x: float(y) for x, y
                        in zip(orbitals_labels[1:], values[1:])}
        return orbitals_ions, orbitals_tot

    @staticmethod
    def get_phase_part(flow_lines, num_atoms):
        """
        軌道の位相の部分を読む
        """
        orbitals_labels = flow_lines.readline().split()
        orbitals_phase = {x: {'r': [], 'i': []} for x in orbitals_labels}
        for _ in range(0, num_atoms):
            real = flow_lines.readline().split()
            imaginaly = flow_lines.readline().split()
            for i in range(0, len(orbitals_labels)):
                orbitals_phase[orbitals_labels[i]]['r'].append(real[i])
                orbitals_phase[orbitals_labels[i]]['i'].append(imaginaly[i])
        return {"orbitals_phase": orbitals_phase}

    def get_data(self, flow_lines):
        """get_1kpointを繰り返す"""
        data = []
        for _ in range(0, self.num_kpoints*2):
            points = self.get_1kpoint(flow_lines,
                                      self.num_bands, self.num_atoms)
            data += points
        return data

    def set_spin_direction(self):
        """
        スピンの向き(up/down)をデータに追記する
        データを半分に割って、上半分にup、残りにdownを振る
        """
        i = len(self.data)
        for i in range(0, int(len(self.data)/2)):
            self.data[i].update({"spin": "up"})
        for i in range(int(len(self.data)/2), len(self.data)):
            self.data[i].update({"spin": "down"})


class ProcarNonC(Procar):
    """
    non-colliner計算のProcarを取り扱う
    """

    def __init__(self, procar):
        Procar.__init__(self, procar)

    def get_weight_part(self, flow_lines, num_atoms):
        """
        軌道のweightの部分を読む
        PROCARは恐らく上から順に
        tot, x, y, z方向のスピン
        """
        meta_orbital_header = re.compile(r"ion\s+s.+tot$")
        line = flow_lines.readline()
        while not meta_orbital_header.match(line):
            line = flow_lines.readline()  # orbital headerまでpass
        orbitals_labels = line.split()
        # 全方向の和
        orb_all = self.get_single_weight(flow_lines, num_atoms,
                                         orbitals_labels)
        # 各方向
        orb_x = self.get_single_weight(flow_lines, num_atoms, orbitals_labels)
        orb_y = self.get_single_weight(flow_lines, num_atoms, orbitals_labels)
        orb_z = self.get_single_weight(flow_lines, num_atoms, orbitals_labels)
        return {"orbitals_ions": orb_all[0], "orbitals_tot": orb_all[1],
                "orbitalx_ions": orb_x[0], "orbitalx_tot": orb_x[1],
                "orbitaly_ions": orb_y[0], "orbitaly_tot": orb_y[1],
                "orbitalz_ions": orb_z[0], "orbitalz_tot": orb_z[1]}

    def get_data(self, flow_lines):
        """get_1kpointを繰り返す"""
        data = []
        for _ in range(0, self.num_kpoints):
            points = self.get_1kpoint(flow_lines,
                                      self.num_bands, self.num_atoms)
            data += points
        return data

    def set_spin_direction(self):
        """
        スピンの向き(up/down)をデータに追記する
        """
        directions = ['orbitalx_tot', 'orbitaly_tot', 'orbitalz_tot']
        weights = [DataBox(self[x]) for x in directions]
        net_weights = (weights[0]['tot'] + weights[1]['tot'] +
                       weights[2]['tot'])
        up_or_down = {True: 'up', False: 'down'}
        spin = [up_or_down[float(x) > 0] for x in net_weights]
        self['spin'] = spin

    def get_sum_energies(self):
        """
        各k点でoccupancyの値を掛けてsumを取る >> PROCAR_bandのoccはずれるので
        Ef以下の場合を積算する
        一方のspin方向のみにtrimしておかないと両方混ざった結果を与える
        """
        results = []
        for i in range(0, self.num_kpoints):
            one_point = DataBox(self.trim_data(self['kpoint_id'], i+1))
            # sum_energies = 0
            # for energy in one_point.data:
            #     if energy['energy'] <= 0:
            #         sum_energies += energy['energy']
            orbitals = DataBox(one_point['orbitals_tot'])
            judge = one_point['energy'] < 0
            # sum_energies = sum(one_point['energy'] * one_point['occupancy']
            #                    * (orbitals['dx2'] + orbitals['dxy']))
            sum_energies = sum(one_point['energy'] * judge
                               * (orbitals['dx2'] + orbitals['dxy']))

            results.append({'kpoint_id': i+1, 'Sum_E': sum_energies})
        return DataBox(results)


class Doscar(DataBox):
    """
    DOSCARを読む
    vaspyに移動する？
    """
    labels_orbital = ["s", "s*", "py", "py*", "pz", "pz*", "px", "px*",
                      "dxy", "dxy*", "dyz", "dyz*", "dz2", "dz2*",
                      "dxz", "dxz*", "dx2", "dx2*"]

    def __init__(self, dos):
        DataBox.__init__(self, [])
        self.dos_lines = Cabinet.read_file(dos)
        self.num_atoms = self.get_num_atoms(self.dos_lines[0])
        (self.num_energy,
         self.fermi_energy) = self.get_parameters(self.dos_lines[5])
        print(self.get_parameters(self.dos_lines[5]))
        self.dos_data = self.__prep_dos_data(self.num_atoms,
                                             self.labels_orbital)

    def get_upspin(self):
        """
        up spinのデータを出力する
        """
        self.output_keys = ["s", "py", "pz", "px", "dxy", "dyz", "dz2",
                            "dxz", "dx2"]

    @staticmethod
    def __prep_dos_data(num_atoms, label_orbital):
        """
        不要?? dos_dataの容器を作成
        """
        dos_data = []
        for i in range(0, num_atoms):
            dos_data.append({})
            dos_data[i].update({"Energy": []})
            for j in range(0, 18):
                dos_data[i].update({label_orbital[j]: []})
        return dos_data

    @staticmethod
    def get_num_atoms(line):
        """
        系の元素数をDOSCAR(の一行目)から読み取る
        """
        part = r"\s*([^\s]+)"
        keywords = part * 4 + r"\s*"
        meta = re.compile(keywords)
        match_line = meta.match(line)
        num_atoms = int(match_line.group(1))
        return num_atoms

    @staticmethod
    def get_parameters(line):
        """
        系のパラメータを読み取る
        エネルギーとフェルミエネルギー
        """
        part = r"\s*([^\s]+)"
        keywords = part * 5 + r"\s*"
        meta = re.compile(keywords)
        match_line = meta.match(line)
        num_energies = int(match_line.group(3))
        fermi_energy = float(match_line.group(4))
        return num_energies, fermi_energy

    def get_dos(self, line, labels_orbital, fermi_energy):
        """
        dosを読む
        """
        part = r"(\s+[^\s]+)"
        keywords = part * 19 + r"\s+"
        meta = re.compile(keywords)
        match_line = meta.match(line)
        if not match_line:
            return
        data = {}
        data.update({'energy': float(match_line.group(1)) - fermi_energy})
        for i in range(2, 20):
            data.update({labels_orbital[i-2]: float(match_line.group(i))})
        self.data.append(data)

    @staticmethod
    def get_maxdos(dos_data, outlabels):
        """
        dosの最大値を読み込む
        規格化などに利用
        """
        maximum = []
        for i in range(0, len(outlabels)):
            maximum.append(max(dos_data[outlabels[i]]))
        return max(maximum)

    @staticmethod
    def write_data(dos_data, num_energies, outlabels, fileout):
        """
        データを出力する
        """
        out = open(fileout, 'w')
        out.write("# ")
        for label in outlabels:
            out.write("%s\t" % (label))
        out.write("\n")
        for i in range(0, num_energies):
            for label in outlabels:
                out.write("%f\t" % (dos_data[label][i]))
            out.write("\n")

    def get_data(self):
        """
        dosのデータを読み取る
        upspinとdownspinに別ける??
        """
        for line in self.dos_lines:
            self.get_dos(line, self.labels_orbital, self.fermi_energy)
        # sumDos = {}
        # sumDos.update({"Energy": []})
        # for i in range(0, 18):
        #     sumDos.update({labelOrbit[i-2]: [0] * num_energies})
        # sumDos["Energy"] = dosData[0]["Energy"]
        # for i in range(0, num_atoms):
        #     for j in range(0, 18):
        #         for k in range(0, num_energies):
        #             sumDos[labelOrbit[j]][k] += dosData[i][labelOrbit[j]][k]
        # def makeDOSdata(plotlistOrbital):
        #     outLabelsDown = ["Energy"]
        #     for list in plotlistOrbital:
        #         outLabelsDown.append("%s*" % (list))
        #     writeData(sumDos, outLabelsDown, "dos_down.dat")
        #     outLabelsUp = ["Energy"]
        #     for list in plotlistOrbital:
        #         outLabelsUp.append("%s" % (list))
        #     writeData(sumDos, outLabelsUp, "dos_up.dat")
        # maximum = getMaxDos(sumDos, labelOrbit)
        # print("Ef =%f" % (fermiEnergy))


class Energy(DataBox):
    """
    EnergyPreと入れ替え
    組成も収容できる
    mae計算を追加したらEnergyと交換
    arrayだとlistデータを扱えなくなるためlist版に変更
    """
    def __init__(self, path_list, poscar='POSCAR', output='OSZICAR'):
        """
        initialize
        """
        if not path_list:
            print("path_list is empty !!!")
        self.filetype = {'poscar': poscar, 'output': output}
        self.path_list = path_list
        self.error_path = []
        DataBox.__init__(self, self.get_data())
        #self.data, self.comp = self.get_data()
        #self.trim_bool(self.data, self.data['bool'], 'energy', 'mag')
        self.output_keys = []

    def set_mae(self, osz1, osz2):
        """
        MAEの計算値をself.data_arrayに追加
        """
        initial = self.filetype
        self.filetype['output'] = osz1
        soc1 = DataBox(self.get_data())
        self.filetype['output'] = osz2
        soc2 = DataBox(self.get_data())
        mae = (soc2['energy'] - soc1['energy']) * (10 ** 6)
        # bools = [soc2['bool'][x] and soc1['bool'][x]
        #          for x in range(0, len(soc1['bool']))]
        self.filetype = initial
        self['mae'] = mae

    def get_data(self):
        """
        POSCAR/CONTCARからc/aとvolume,
        OSZICAR/OUTCARからenergyとmagを読む
        元素種、数も読めるようにしたい
        各パラメータのdictをreturn
        """
        data_box = []
        for dirc in self.path_list:
            path_oszi = os.path.join(dirc, self.filetype['output'])
            path_posc = os.path.join(dirc, self.filetype['poscar'])
            result = self.get_data_single(path_posc, path_oszi)
            if result:
                data_box.append(result)
        #data_box.sort(key=lambda x: x['elements'][2])
        return data_box

    def get_data_single(self, poscar='POSCAR', oszicar='OSZICAR'):
        """
        src_dir中のデータを修得
        POSCARからc/aとvolumeとformula
        OSZICAR/OUTCARから読み取ったenergyとmagをreturn
        volume, energy, magの3つは単位原子あたりに換算
        dict形式で出力
        """
        path = os.path.dirname(poscar)
        cova, volume, sum_atoms, elements, num_atoms, spg, spg_num = \
            self.read_posc(poscar)
        energy, mag, valid = self.get_enemag(oszicar)
        try:
            result = {'c/a': cova, 'energy': energy / sum_atoms,
                      'volume': volume / sum_atoms, 'mag': mag / sum_atoms,
                      'elements': elements, 'num_atoms': num_atoms,
                      'bool': valid, 'path': path,
                      'spg': spg, 'spg_num': spg_num}
            return result
        except TypeError:
            print(cova, volume, sum_atoms, elements, num_atoms, energy, mag,
                  valid)

    def read_posc(self, poscar='POSCAR'):
        """
        Read cova, volume, num_atoms from POSCAR file.
        spg も読めるように変更
        """
        posc = vaspy.Poscar(poscar)
        latt_a, _, latt_c = posc.get_lattice_length()[0:3]
        cova = latt_c/latt_a
        volume = posc.get_cell_volume()
        elements = posc.elements
        num_atoms = posc.num_atoms
        sum_atoms = sum(num_atoms)
        spg, spg_num = self.print_spg(poscar)
        return cova, volume, sum_atoms, elements, num_atoms, spg, spg_num

    @staticmethod
    def print_spg(src='POSCAR'):
        """
        space group を return
        """
        srcpos = Poscar.from_file(src)
        finder = SpacegroupAnalyzer(srcpos.structure,
                                    symprec=5e-2, angle_tolerance=8)
        spg = finder.get_spacegroup_symbol()
        spg_num = finder.get_spacegroup_number()
        return spg, spg_num

    @staticmethod
    def get_enemag(oszicar='OSZICAR'):
        """
        Read energy and magnetic momentum from OSZICAR/OUTCAR.
        Judge file type from first line of the read file.
        """
        try:
            lines = Cabinet.read_file(oszicar)
        except IOError:
            lines = ['error']
        try:
            head = lines[0].split()[0]
        except IndexError:
            print(oszicar)
            exit()
        if head == 'N':  # OSZICAR
            output = vaspy.Oszicar(oszicar)
        elif head[0:4] == 'vasp':  # OUTCAR
            output = vaspy.Outcar(oszicar)
        else:
            print("I cannot read {0}".format(oszicar))
            return None, None, False
        if not output.results:
            print(oszicar)
            return None, None, False
        energy = output.results[-1]['energy']
        mag = output.results[-1]['mag']
        return energy, mag, True

    def __getitem__(self, key):
        array = [x[key] for x in self.data]
        return np.array(array) #pylint: disable=E1101

    def __setitem__(self, key, array):
        for i in range(0, len(self.data)):
            if len(self.data) != len(array):
                print("Dimension of data is different")
                break
            try:
                self.data[i][key] = array[i]
            except KeyError:
                self.data[i].update({key: array[i]})

    def set_enthalpy(self):
        """
        Convert energy (eV/atom) to enthalpy (kJ/mol)
        [solid.REFERENCE[x] for x in data['elements'][:,0]]
        のような記述も出来るけどとりあえずは地味に
        """
        for data in self.data:
            elements = data['elements']
            num_atoms = data['num_atoms']
            sum_atoms = sum(num_atoms)
            ref_energy = 0
            for element, num in zip(elements, num_atoms):
                ref_energy += solid.REFERENCE[element] * num
            ref_energy = ref_energy / sum_atoms
            enthalpy = ((data['energy'] - ref_energy) * 96.485344520851)
            data.update({'enthalpy': enthalpy})

    def set_elements_z(self):
        """
        原子番号に読みかえる
        value_z = [[solid.ELEMENTS[x]['Z'] for x in y]
                    for y in self['elements']]
        """
        for data in self.data:
            value_z = [solid.ELEMENTS[x]['Z'] for x in data['elements']]
            data.update({'Z': value_z})

    def set_comp_dict_i(self):
        """
        組成比を整数で作成
        """
        comp_all = set([x for y in self['elements'] for x in y])
        for data in self.data:
            comp_dict = {elem: 0 for elem in comp_all}
            comp_dict_r = {elem: num for elem, num
                           in zip(data['elements'], data['num_atoms'])}
            comp_dict.update(comp_dict_r)
            data.update({'comp_dict_i': comp_dict})

    def set_comp_dict_f(self):
        """
        組成比を分率で作成
        """
        comp_all = set([x for y in self['elements'] for x in y])
        for data in self.data:
            comp_dict = {elem: 0 for elem in comp_all}
            sum_atoms = float(sum(data['num_atoms']))
            comp_dict_r = {elem: num/sum_atoms for elem, num
                           in zip(data['elements'], data['num_atoms'])}
            comp_dict.update(comp_dict_r)
            data.update({'comp_dict_f': comp_dict})

    def print_comp(self, *elements):
        """
        convex hull用
        """
        for data in self.data:
            sum_atoms = sum(data['num_atoms'])
            comp = []
            for element in elements:
                try:
                    i = data['comp_dict'][element]
                except KeyError:
                    i = 0
                comp.append(str(i/sum_atoms))
            comp.append(str(data['enthalpy']))
            line = " ".join(comp)
            print(line)

    def __str__(self):
        """
        Print self.data
        特定の要素のみを書き出したい場合は
        self.output_keysに指定する
        defaultでは全てを書き出す
        """
        keys = self.output_keys
        if not keys:
            keys = self.data[0].keys()
        types_dict = [x for x in keys if type(self.data[0][x]) is dict]
        out_lines = "\t".join(keys) + "\n"
        for key in types_dict:
            head = "".join(key.split('_dict')) + "_"
            label = [head + str(x) for x in sorted(self.data[0][key])]
            alt_label = "\t".join(label)
            out_lines = out_lines.replace(key, alt_label)
        for data in self.data:
            line = []
            for key in keys:
                if type(data[key]) is dict:
                    data[key] = "\t".join([str(data[key][x]) for x
                                           in sorted(data[key])])
                line.append(str(data[key]))
            out_lines += "\t".join(line) + "\n"
        return out_lines

    def separate_data(self, label, pos):
        """
        self.data の self.data[label] が同一の要素毎に分割する
        分割の判定要素が list[pos] のケース
        plot 用
        """
        tmp_values = [x[pos] for x in self[label]]
        values = sorted(set(tmp_values), key=tmp_values.index)
        separated_data = []
        for val in values:
            data = [x for x in self.data if x[label][pos] == val]
            separated_data.append(DataBox(data))
        return separated_data

    def table_combi(self, sites, out_key):
        """
        self.dataをcombi用のtableに整形
        POSCARでのcombiの元素位置をsiteに指定
        __str__()でx:site[0], y:site[1]の元素でtableが出力
        """
        data_list = self.separate_data('elements', sites[1])
        tmp_table = []

        tmp_key = [x[sites[0]] for x in self['elements']]
        keys = sorted(set(tmp_key), key=tmp_key.index)

        for data in data_list:
            add_dict = {out_key: data.data[0]['elements'][sites[1]]}
            add_dict.update({x['elements'][sites[0]]: x[out_key]
                             for x in data.data})
            tmp_table.append(add_dict)
        table = DataBox(tmp_table)
        table.output_keys = [out_key] + keys
        print(table.output_keys)
        return table

    def set_formula(self):
        """
        pylabで表示するformulaをsetする
        """
        for data in self.data:
            formula = ""
            for elem, i in zip(data['elements'], data['num_atoms']):
                formula += "{0}$_{1}$".format(elem, i)
            data.update({'formula': formula})

    def alt_elements_from_dir(self):
        """
        組成式の入ったelem_...形式のディレクトリ名から
        elements & num_atoms を再読み込みする
        """
        for data in self.data:
            formula = os.path.basename(data['path']).split('_')[-1]
            upper = len(re.split('[A-Z]', formula)) - 1
            meta = r'([A-Z][a-z]*)(\d*)' * upper
            match = re.match(meta, formula)
            elements = [match.group(2*i+1) for i in range(upper)]
            num_atoms = []
            for i in range(upper):
                try:
                    num_atoms.append(int(match.group(i*2+2)))
                except ValueError:
                    num_atoms.append(1)
            ratio = sum(num_atoms) / sum(data['num_atoms'])
            norm_num_atoms = [x / ratio for x in num_atoms]
            data['elements'] = elements
            data['num_atoms'] = norm_num_atoms

    def alt_cova_posstd(self):
        """
        poscat_std から c/a を再読みこみする
        """
        for data in self.data:
            path = os.path.join(data['path'], 'POSCAR_std')
            posc = vaspy.Poscar(path)
            latt_a, _, latt_c = posc.get_lattice_length()[0:3]
            cova = latt_c/latt_a
            data['c/a'] = cova

    def set_mag_each_site(self):
        """
        各サイトごとの磁気モーメントを習得
        label を設定する
        """
        for data in self.data:
            path = os.path.join(data['path'], 'POSCAR')
            with open(path, 'r') as rfile:
                lines = rfile.readlines()
            elements = lines[5].split()
            nums = lines[6].split()
            total = sum([int(x) for x in nums])
            symbols = [i+str(x+1) for i, j in zip(elements, nums)
                       for x in range(int(j))]
            path = os.path.join(data['path'], 'OUTCAR')
            with open(path, 'r') as rfile:
                lines = rfile.readlines()
            p = lines.index(" magnetization (x)\n")
            mag = [float(lines[i].split()[4]) for i in range(p+4, p+4+total)]
            data['symbol_mag/site'] = symbols
            data['mag/site'] = mag

class MurnaghanFit(Energy):
    """
    voldepを計算した結果をMurnaghan_fitして収集する
    """
    def __init__(self, dir_list):
        Energy.__init__(self, dir_list)

    def get_data(self):
        """
        fitting resultからE0, B, B1, volume, それらのfitting誤差を読む
        """
        data_box = []
        for path in self.path_list:
            print(path)
            fit_res = self.fit_murnaghan(path)
            if not fit_res:
                continue
            path_posc = os.path.join(path, self.filetype['poscar'])
            cova, _, _, elements, num_atoms = self.read_posc(path_posc)
            data_box.append({'c/a': cova, 'elements': elements,
                             'num_atoms': num_atoms, 'path': path,
                             'energy': fit_res[0][0], 'errE': fit_res[1][0],
                             'B': fit_res[0][1], 'errB': fit_res[1][1],
                             'B1': fit_res[0][2], 'errB1': fit_res[1][2],
                             'volume': fit_res[0][3], 'errV': fit_res[1][3]})
        #data_box.sort(key=lambda x: x['elements'][2])
        return data_box

    def fit_murnaghan(self, path):
        osz_path = os.path.join(path, 'voldep/volume_*/OSZICAR')
        dir_list = [os.path.dirname(x) for x in glob.glob(osz_path)]
        if not dir_list:
            self.error_path.append(path)
            return None
        data = Energy(dir_list, 'POSCAR', 'OSZICAR')
        data.data.sort(key=lambda x: x['volume'])
        if len(data.data) <= 4:  # fitting paramよりデータが少ない
            self.error_path.append(path)
            return None
        fit_res = FitData.Murnaghan_fit(data['volume'], data['energy'])
        if fit_res[3].all() == 0:
            self.error_path.append(path)
            return None
        lines = ""
        for i in range(0, len(fit_res[0])):
            lines += "{0}\t{1}\n".format(fit_res[0][i], fit_res[1][i])
        fname = os.path.join(path, 'fit_curve.dat')
        with open(fname, 'w') as wfile:
            wfile.write(lines)
        pylab.plot(data['volume'], data['energy'])
        pylab.plot(fit_res[0], fit_res[1])
        fname = os.path.join(path, 'murnaghan_plot.eps')
        pylab.savefig(fname)
        pylab.close()
        lines = ""
        lines += "Energy\t{0[2][0]}\terrE\t{0[3][0]}\n".format(fit_res)
        lines += "B\t{0[2][1]}\terrB\t{0[3][1]}\n".format(fit_res)
        lines += "B1\t{0[2][2]}\terrB1\t{0[3][2]}\n".format(fit_res)
        lines += "Volume\t{0[2][3]}\terrV\t{0[3][3]}\n".format(fit_res)
        fname = os.path.join(path, 'fit_results.dat')
        with open(fname, 'w') as wfile:
            wfile.write(lines)
        return fit_res[2:]



class EnergyPre(object):
    """
    energy関係のoutputを収集する
    pathのlistを引数として、
    list中の全てのposcar, oszicarファイルを読む
    c/a, volume, energy, magのself.data (array)を作成する
    volume, energy, magは単位原子あたりに換算
    データ形式をnp.arrayではなくリスト形式に変更する
    EnergyList(object)

    valid=0/1の項も追加
    valid=1のデータにはnoneを与える

    ↑　必要か??
    各データの次元を揃える為にあった方がいい
    使う方向で考える
    """
    def __init__(self, path_list, poscar='POSCAR', output='OSZICAR'):
        """
        initialize
        """
        self.filetype = {'poscar': poscar, 'output': output}
        self.path_list = path_list
        self.data = self.get_data()
        #self.data, self.comp = self.get_data()
        self.trim_bool(self.data, self.data['bool'], 'energy', 'mag')
        self.output_list = ['c/a', 'energy', 'mag', 'volume']

    def alt_filetype(self, poscar=None, output=None):
        """
        filetypeを変更
        """
        if poscar:
            self.filetype['poscar'] = poscar
        if output:
            self.filetype['output'] = output

    @staticmethod
    def trim_bool(array, bools, *labels):
        """
        boolsがFalseの場合array[label]にNoneを代入
        """
        for i in range(0, len(bools)):
            if not bools[i]:
                for label in labels:
                    array[label][i] = None

    def get_data(self):
        """
        POSCAR/CONTCARからc/aとvolume,
        OSZICAR/OUTCARからenergyとmagを読む
        元素種、数も読めるようにしたい
        indexをkeyにvalueにnp.arrayのdictをreturn
        """
        data_box = []
        for dirc in self.path_list:
            print(dirc)
            path_oszi = os.path.join(dirc, self.filetype['output'])
            path_posc = os.path.join(dirc, self.filetype['poscar'])
            result = self.get_data_single(path_posc, path_oszi)
            data_box.append(result)
        data_box.sort(key=lambda x: x[0])
        dtype = [('c/a', 'f8'), ('energy', 'f8'),
                 ('volume', 'f8'), ('mag', 'f8'), ('bool', 'bool')]
        return np.array([tuple(x) for x in data_box], dtype=dtype) #pylint: disable=E1101

    def get_data_single(self, poscar='POSCAR', oszicar='OSZICAR'):
        """
        src_dir中のデータを修得
        POSCARからc/aとvolumeとformula
        OSZICAR/OUTCARから読み取ったenergyとmagをreturn
        volume, energy, magの3つは単位原子あたりに換算
        """
        cova, volume, num_atoms = self.get_lattices(poscar)
        energy, mag, valid = self.get_enemag(oszicar)
        result = ([cova, energy / num_atoms,
                   volume / num_atoms, mag / num_atoms,
                   valid])
        return result

    @staticmethod
    def get_lattices(poscar='POSCAR'):
        """
        Read cova, volume, num_atoms from POSCAR file.
        """
        posc = vaspy.Poscar(poscar)
        latt_a, _, latt_c = posc.get_lattice_length()[0:3]
        cova = latt_c/latt_a
        volume = posc.get_cell_volume()
        num_atoms = sum(posc.num_atoms)
        return cova, volume, num_atoms

    @staticmethod
    def get_enemag(oszicar='OSZICAR'):
        """
        Read energy and magnetic momentum from OSZICAR/OUTCAR.
        """
        try:
            lines = Cabinet.read_file(oszicar)
        except IOError:
            lines = ['error']

        if lines[0].split()[0] == 'N':  # OSZICAR
            output = vaspy.Oszicar(oszicar)
        elif lines[0].split()[0][0:4] == 'vasp':  # OUTCAR
            output = vaspy.Outcar(oszicar)
        else:
            print("I cannot read {0}".format(oszicar))
            return 0, 0, False
        if not output.results:
            print(oszicar)
            return 0, 0, False
        energy = output.results[-1]['energy']
        mag = output.results[-1]['mag']
        return energy, mag, True

    def get_mae(self, osz1, osz2, label):
        """
        MAEの計算値をself.data_arrayに追加
        """
        initial = self.filetype
        self.alt_filetype(output=osz1)
        soc1 = self.get_data()
        self.alt_filetype(output=osz2)
        soc2 = self.get_data()
        mae = (soc2['energy'] - soc1['energy']) * (10 ** 6)
        bools = [soc2['bool'][x] and soc1['bool'][x]
                 for x in range(0, len(soc1['bool']))]
        dtype = (label, 'f8')
        self.data = Array.add_data(self.data, mae, dtype)
        self.trim_bool(self.data, bools, 'mae')
        self.filetype = initial
        line = "MAE = '{0}' - '{1}'".format(osz2, osz1)
        print(line)
        self.output_list.append('mae')

    @staticmethod
    def get_composition(poscar, potcar):
        """
        POSCARとPOTCARから元素名と構成数を読み込み
        vaspのバージョンが5の場合POSCARのみから読む
        """
        pos_obj = vaspy.Poscar(poscar)
        if pos_obj.vasp_version == 4:
            elements = vaspy.Potcar.get_composition(potcar)
        else:
            elements = pos_obj.elements
        num_atoms = pos_obj.num_atoms
        return elements, num_atoms

    def __str__(self):
        """Print data in self.output_list"""
        out_lines = "\t".join(self.output_list) + "\n"
        for i in range(0, self.data.shape[0]):
            line = []
            for output in self.output_list:
                line.append(str(self.data[output][i]))
            out_lines += "\t".join(line) + "\n"
        return out_lines

    def get_enthalpy(self, elements, num_atoms):
        """Convert energy (eV/atom) to enthalpy (kJ/mol)"""
        ref_energy = 0
        total_num = sum(num_atoms)
        for element, num in zip(elements, num_atoms):
            ref_energy += solid.REFERENCE[element] * num
        ref_energy = ref_energy / total_num
        enthalpy = ((self.data['energy'] - ref_energy) * 96.485344520851)
        dtype = ('enthalpy', 'f8')
        self.data = Array.add_data(self.data, enthalpy, dtype)
        self.output_list.append('enthalpy')
