#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
cvm imput file を修正したりする
"""
import re
from fractions import Fraction


class CVMEnergies(object):
    """
    CVMEnergy の集合を list で管理
    """
    def __init__(self, energies):
        self.energies = energies

    @classmethod
    def from_energies_file(cls, path):
        """
        energies.txt から object を生成
        """
        with open(path, 'r') as rfile:
            lines = rfile.readlines()
        energies = [CVMEnergy.from_str(x) for x in lines]
        return cls(energies)

    def __str__(self):
        return "\n".join([str(x) for x in self]) + "\n"

    def __getitem__(self, x):
        return self.energies[x]

    def exchange_elements(self, elem1, elem2):
        """
        元素 elem1 と elem2 とを入れ替えた energy label を作成する
        """
        for energy in self:
            try:
                energy.label.exchange_elements(elem1, elem2)
            except ValueError:
                line = str(energy.label)
                line += " might be disorder str >>> It has been skipped !!!"
                print(line)

    def make_file(self, path):
        """
        energies.txt を作成する
        """
        with open(path, 'w') as wfile:
            wfile.write(str(self))


class CVMEnergy(object):
    """
    CVM 形式の energy
    """
    def __init__(self, label, energy, unit):
        self.label = label
        self.energy = energy
        self.unit = unit

    def __str__(self):
        return str(self.label) + " 1 " + str(self.energy) + " " + str(self.unit)

    @classmethod
    def from_str(cls, strings):
        """
        str から object を生成
        """
        label = CVMLabel.from_str(strings.split()[0])
        energy = float(strings.split()[2])
        unit = int(strings.split()[3])
        return cls(label, energy, unit)

    def exchange_elements(self, elem1, elem2):
        """
        元素のラベルを入れ替える
        """
        self.label.exchange_elements(elem1, elem2)


class CVMStrs(object):
    """
    CVMStructure の集合を list で管理
    """
    def __init__(self, structures):
        self.strs = structures

    @classmethod
    def from_str_file(cls, path):
        """
        **.str ファイルから object を生成
        """
        with open(path, 'r') as rfile:
            lines = rfile.readlines()
        meta = re.compile(r".*NPH=.*")
        lines.append(lines[0])
        i = 0
        strs = []
        while i < len(lines) - 1:
            j = i
            i += 1
            while not meta.match(lines[i]):
                i += 1
            strs.append(CVMStructure.from_lines(lines[j:i]))
        return cls(strs)

    def __getitem__(self, x):
        return self.strs[x]

    def __str__(self):
        return "".join([str(x) for x in self])

    def exchange_elements(self, elem1, elem2):
        """
        元素 elem1 と elem2 とを入れ替えた structure を作成する
        """
        for structure in self:
            try:
                structure.exchange_elements(elem1, elem2)
            except ValueError:
                line = str(structure.label)
                line += " might be disorder str >>> It has been skipped !!!"
                print(line)

    def make_file(self, path):
        """
        **.str を作成する
        """
        with open(path, 'w') as wfile:
            wfile.write(str(self))

    def add_site(self, primitive_z, origin, trans, label, insert_idx):
        """
        site を追加する
        """
        for structure in self:
            structure.add_site(primitive_z, origin, trans, label, insert_idx)


class CVMStructure(object):
    """
    CVM 形式の structure
    """
    def __init__(self, label, brav, disord, nph, lattice, sites):
    #pylint: disable=R0913
        self.label = label
        self.brav = brav
        self.disord = disord
        self.nph = nph
        self.lattice = lattice
        self.sites = sites

    @classmethod
    def from_lines(cls, lines):
        """
        lines から object を return
        """
        label = CVMLabel.from_str(lines[0].split()[0])
        brav = lines[0].split()[1].split("=")[-1]
        disord = lines[0].split()[2].split("=")[-1]
        nph = int(lines[0].split()[3].split("=")[-1])

        if len(lines[1].split("=")) <= 2:
            lattice = CVMLattice(
                [[Fraction(x) for x in line.split("=")[-1].split()]
                 for line in lines[1:4]])
        else:
            lattice = CVMLattice([[Fraction(x) for x in y.split()[0:3]]
                                  for y in lines[1].split("=")[1:]])
        sites = []
        for line in lines[4:]:
            atom = line.split("=")[-1].split()[0]
            coord = [Fraction(x) for x in line.split("=")[-1].split()[1:]]
            sites.append(CVMSite(atom, coord))
        args = {'label': label, 'brav': brav, 'disord': disord,
                'nph': nph, 'lattice': lattice, 'sites': sites}
        return cls(**args) #pylint: disable=W0142

    def add_site(self, primitive_z, origin, trans, label, insert_idx):
    #pylint: disable=R0913
        """
        サイトを追加する
        複製元 (origin) に trans vector を足してサイトを追加する
        primitive cell 中の原子数を指定することで、
        supercell の等価なサイトに対して同じ処理を繰り返す
        args:
            primitive_z: primitive cell 中の原子数 (int)
            origin: 複製元の site index (int)
            trans: trans vector (3d-list)
            label: 追加するサイトの atom label
        """
        size = len(self.sites) / primitive_z
        for i in range(size):
            new_coord = [x+y for x, y in
                         zip(self.sites[origin+i*primitive_z+i].coord, trans)]
            new_site = CVMSite(label, new_coord)
            self.sites.insert(insert_idx+i*primitive_z+i, new_site)
        idx = self.label.elements.index(label[0])
        self.label.num_atoms[idx] += size

    def __str__(self):
        lines = ""
        lines += str(self.label) + " "
        lines += "BRAV=" + self.brav + " "
        lines += "DISORD=" + self.disord + " "
        lines += "NPH=" + str(self.nph) + "\n"
        lines += str(self.lattice)
        for i in range(len(self.sites)):
            lines += "   " + "ATOM" + str(i+1) + "="
            lines += str(self.sites[i])
        return lines

    def exchange_elements(self, elem1, elem2):
        """
        elem1, elem2 にマッチする原子 label を交換する
        """
        newlab = {elem1: elem2, elem2: elem1}
        for site in self.sites:
            if site.atom in [elem1, elem2]:
                site.atom = newlab[site.atom]
        self.label.exchange_elements(elem1, elem2)

    def get_conc(self, elem):
        """
        elem の濃度を return する
        """
        return self.label.get_conc(elem)


class CVMLabel(object):
    """
    CVM 形式の構造 label
    """
    def __init__(self, label, elements, num_atoms):
        self.label = label
        self.elements = elements
        self.num_atoms = num_atoms
        self._formula = self.formula

    def get_conc(self, elem):
        """ elem の濃度を return する"""
        idx = self.elements.index(elem)
        return self.num_atoms(idx) / sum(self.num_atoms)

    @property
    def formula(self):
        """
        組成式を return
        """
        self._formula = "".join(
            ["".join([x, str(y)])
             for x, y in zip(self.elements, self.num_atoms)])
        return self._formula

    @classmethod
    def from_str(cls, line):
        """
        line から object を作成
        """
        label = line.split()[0].split("*")[0]
        formula = line.split()[0].split("*")[1]
        elements = [x for x in formula if x in "ABCDEFGHIJKLMNOPQRSTUVWXY"]
        num_atoms = []
        for i in range(len(elements)):
            head = formula.index(elements[i])
            tail = (formula + "Z").index((elements+["Z"])[i+1])
            try:
                num_atoms.append(int(formula[head+1:tail]))
            except ValueError:
                num_atoms.append(1)
        return cls(label, elements, num_atoms)

    def exchange_elements(self, elem1, elem2):
        """
        elem1, elem2 に該当する元素数を入れ替える
        """
        at1 = elem1[0]
        at2 = elem2[0]
        pos1 = self.elements.index(at1)
        pos2 = self.elements.index(at2)
        num1 = self.num_atoms[pos1]
        num2 = self.num_atoms[pos2]
        self.num_atoms[pos2] = num1
        self.num_atoms[pos1] = num2

    def __str__(self):
        return self.label + "*" + self.formula


class CVMLattice(object): #pylint: disable=R0903
    """
    CVM 形式の lattice の class
    """
    def __init__(self, lattice):
        self.a = lattice[0]
        self.b = lattice[1]
        self.c = lattice[2]
    def __str__(self):
        lines = ""
        lines += "   " + "A_CART="
        lines += " ".join(str(Fraction(x)) for x in self.a) + "\n"
        lines += "   " + "B_CART="
        lines += " ".join(str(Fraction(x)) for x in self.b) + "\n"
        lines += "   " + "C_CART="
        lines += " ".join(str(Fraction(x)) for x in self.c) + "\n"
        return lines


class CVMSite(object): #pylint: disable=R0903
    """
    CVM 形式の site の class
    """
    def __init__(self, atom, coord):
        self.atom = atom
        self.coord = coord

    def __str__(self):
        lines = ""
        lines += self.atom + " "
        lines += " ".join(str(Fraction(x)) for x in self.coord) + "\n"
        return lines


