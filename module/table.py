#!/opt/local/bin/python
# -*- coding: utf-8 -*-
"""
energy, c/a, mag, volumeなどの情報を読んでdata_arrayを作成する
→ collection.pyとする
"""
import math
import os
import commands
import numpy
import pylab
import scipy.optimize
import itertools
import glob
from makeSeries import MakePattern, Combinatorial
from commopy import Cabinet, Bash
import solid
import vaspy
import grapy

Reference = {"Fe": -8.30986305725147, "Co": -7.10646786208887,
             "Mn": -9.04116583023117, "Ni": -5.57058693907496,
             "Cr": -9.64529903086873, "V": -9.08055004986443,
             "Ti": -7.89042111815724, "Cu": -3.71797397606303,
             "Zn": -1.26775730145801, "Ga": -3.0281282763815,
             "Ge": -4.62163453127987, "In": -3.18887807809515,
             "Sn": -4.00527149401934, "Pd": -5.18089045007217,
             "Rh": -7.26929051746743, "Pt": -6.05335455763213,
             "Ru": -9.20388484292648, "Nb": -10.0945900627376,
             "Zr": -8.54761937261647, "Al": -3.74129357979325,
             "Si": -5.42476168374865, "Sb": -4.12207679803694,
             "Pb": -3.70719513700551, "N": -8.31608161824619,
             "B": -6.67811265223554, "C": -9.1094663584899,
             "F": -1.85748416912446, "P": -5.37439133969137,
             "S": -4.25214836550994, "Cl": -1.78700714451757}

CombiPara = [["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu"],
             ["Al", "Si", "Zn", "Ga", "Ge", "In", "Sn", "Sb", "Pb"]]
CombiPara = [["Ti", "V", "Cr", "Mn", "Co", "Ni", "Cu"],
             ["Al", "Si", "Zn", "Ga", "Ge", "In", "Sn", "Sb", "Pb"]]
#CombiPara = [["Co", "Ni", "Cu"],
#             ["Al", "Si", "Zn", "Ga", "Ge", "In", "Sn", "Sb", "Pb"]]
#CombiPara = [["Co"],["Al","Ga","Ge","Si", "Sn"]]
#CombiPara = [["Cu"],["Al","Ga"]]
CombiPara = [["Cu"], ["Ga"]]
Dir_List = ["fixed_1", "fixed_2", "fixed_3", "fixed_4", "fixed_5", "fixed_6",
            "fixed_7", "fixed_8", "fixed_9", "fixed_10", "fixed_11",
            "fixed_12", "fixed_13", "fixed_14", "fixed_15", "fixed_20",
            "fixed_21", "fixed_22", "fixed_23", "fixed_24", "fixed_25",
            "fixed_26"]


def main():
    #plt_akaikkr()

    #X2FeZ
    os.chdir("/Users/enoki/Documents/01_ResearchData/"
             "Calculations/01_VASP/13_Herusler/fixed_lattice")
    plt_double('typeA_MAE', 'typeB_MAE', '01_output', 8)
    #return

    #  akai_kkr
    #os.chdir("/Users/enoki/Documents/01_ResearchData/"
    #         "Calculations/06_AKAI-KKR/00_check/Compare_Enthalpy/01_AlNi/")
    #akai_voldep("AlA1", "AlA2", "AlNi", "NiA1", "NiA2")
    return

    os.chdir("/Users/enoki/Documents/01_ResearchData/"
             "Calculations/01_VASP/13_Herusler/fixed_lattice")
    #plt_single('typeA_60GPa', '02_output_60GPa', 8)
    plt = PlotMAE(HeuslerFe2YZ)
    plt.plt_double_MAE('typeA_MAE', 'typeB_MAE', '04_output_MAE')
    return
    os.chdir("/Users/enoki/Documents/01_ResearchData/"
             "Calculations/06_AKAI-KKR/01_Fe2YZ_Heuslar/CoAl_gga91/")
    dir_list = commands.getoutput('ls').split()
    data = []
    for dir in dir_list:
        os.chdir(dir)
        a = plt_akaikkr('../../01_output_gga91', dir)
        data.append(a)
        os.chdir('..')
    print data
    #ar = numpy.array([1, 2, 3, 4, 5])
    #ar2 = numpy.array([True, True, False, True, False])
    #print ar
    #print ar[ar2]


def akai_voldep(*dir_list):
    i = 0
    output = []
    for dir in dir_list:
        os.chdir(dir)
        plt = PlotResults(dir, 1)
        Data = AkaiKKR()
        HD = MakePattern([[]], "")
        dirlist = HD.get_keyword_dir(".", "latt*")
        Data.data_array = Data.get_data(dirlist, 'output')
        fitpara = plt.Ene_Vol_single(Data.data_array, i, '')
        out = [dir] + fitpara.tolist()
        output.append(out)
        #print a
        #xi, yi, Murnaghan_param = plt.Murnaghan_fit(Data.data_array['volume'],
        #                                            Data.data_array['energy'])
        #plt.Ene_CovA(Data.data_array, 1, "none")
        #print Murnaghan_param
        i += 5
        os.chdir('..')
        pylab.savefig('%s.eps' % dir)
        #pylab.show()
    for out in output:
        for val in out:
            print val,
        print '\n',


def plt_akaikkr():
    data = AkaiKKR()
    path = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
            "06_AKAI-KKR/01_Fe2YZ_Heuslar/")

    comp_dir_list = glob.glob(os.path.join(path, '*',
                              'Conc_100/latt_7.60/CovA_1.00/output'))
    comp_dir_list = [x.split('/Conc_100')[0] for x in comp_dir_list]
    comp_dir_list = [os.path.join(path, 'CoAl-CoV')]
    conc_dir_list = glob.glob(os.path.join(comp_dir_list[0], '*'))
    latt_dir_list = glob.glob(os.path.join(conc_dir_list[0], '*'))
    comp_dir_list = [os.path.basename(x) for x in comp_dir_list]
    conc_dir_list = [os.path.basename(x) for x in conc_dir_list]
    latt_dir_list = [os.path.basename(x) for x in latt_dir_list]

    def plt_one(comp, conc, plt):
        print comp
        print conc
        dir_list = [os.path.join(path, comp, conc, x, 'CovA_1.00')
                    for x in latt_dir_list]
        data_array = data.get_data(dir_list, 'output')
        i = 0
        print data_array['energy']
        xi, yi, Murnaghan_param = plt.Murnaghan_fit(data_array['volume'],
                                                    data_array['energy'])
        param = Murnaghan_param
        return param

        pltstyle_dict = plt.mkstyle(i)
        pltstyle_dict.update({'linestyle': 'None'})
        data_label = 'Fe$_2$(Co,V)Al'

        x, labelx = data_array['volume'], 'c/a'
        y1, label1 = data_array['energy'], 'Energy (eV/atom)'
        y2, label2 = data_array['mag'], 'Mag. ($\mu_B$/atom)'

        plt.ax1.plot(x, y1, label=data_label, **pltstyle_dict)
        pylab.ylabel(label1, fontsize=20, horizontalalignment='left')
        pylab.yticks(fontsize=13)
        plt.ax1.yaxis.set_label_coords(-.12, 0.5)
        plt.ax1.legend()

        plt.ax2 = pylab.subplot(212)
        plt.ax2.plot(x, y2, **pltstyle_dict)
        pylab.ylabel(label2, fontsize=20)
        pylab.xlabel(labelx, fontsize=20)
        pylab.xticks(fontsize=13)

        style = {}
        style.update({'color': pltstyle_dict['color']})
        plt.add_plt(plt.ax1, xi, yi, style)
        pylab.show()

        i += 3
        return param
    #dtype = [('c/a', 'f8'), ('E0', 'f8')]

    comp = comp_dir_list[0]

    def plt_comp(comp):
        data_conc = []
        data_E0 = []
        data_V0 = []
        print(comp)
        plt = PlotResults('', 1)
        for conc in conc_dir_list:
            E0, _, _, V0 = plt_one(comp, conc, plt)
            conc_val = float(conc.split('_')[-1])
            data_conc.append(conc_val)
            data_E0.append(E0)
            data_V0.append(V0)
        delta = data_E0[data_conc.index(100)] - data_E0[data_conc.index(0)]
        col = [data_E0[i] - data_E0[data_conc.index(0)] -
               delta * data_conc[i] / 100
               for i in range(0, len(data_conc))]
        plt.ax1 = pylab.subplot(111)
        label = comp.split('-')
        label = "Fe$_2${0[0]} - Fe$_2${0[1]}".format(label)
        plt.ax1.plot(data_conc, col, 'o', label=label)
        plt.ax1.legend()
        pylab.ylabel('delta', fontsize=20)
        pylab.xlabel('concentration', fontsize=20)
        dst = ("/Users/enoki/Documents/01_ResearchData/Calculations/"
               "06_AKAI-KKR/02_antiperovskite/01_Nitrogen_results/%s" % comp)
        pylab.savefig('%s.eps' % dst)
    print(len(comp_dir_list))
    i = 0
    for comp in comp_dir_list[i:]:
        print(i)
        plt_comp(comp)
        i += 1


def plt_single_old(in_dir, out_dir, atom_num):
    """
    2014/01/28 energy, mom, volumeのc/a依存性をプロット
    """
    Combi = MakePattern(CombiPara, '')
    print Combi.compos
    for compo in Combi.compos:
        print compo
        dir_list = [os.path.basename(x) for x in
                    glob.glob(os.path.join(in_dir, compo, 'fixed_*'))]
        #Combi.mkDirList(in_dir + '/' + compo + '/', dir_list)
        all_dir_list = [os.path.join(in_dir, compo, x) for x in dir_list]
        DataA = HeuslerFe2YZ(all_dir_list, 'POSCAR', 'OSZICAR',
                             'regular')
        #DataA.get_enthalpy(*Combi.compos[compo])

        plt = PlotResults('Fe$_2$%s' % compo)
        plt.EneMomVol_CovA(DataA.data_array, 'blue', DataA.label)
        pylab.savefig('%s/%s.eps' % (out_dir, compo))


class PlotMAE(object):
    def __init__(self, data_object):
        #Combi = MakePattern(CombiPara, '')
        self.compos = ['CoAl', 'CoGa', 'CoGe', 'CoSi', 'CoSn', 'CuAl', 'CuGa',
                       'CuGe', 'CuSi', 'FeAl', 'FeGa', 'FeGe', 'FeSi', 'NiAl',
                       'NiGa', 'NiGe', 'NiSi', 'NiSn', 'NiZn', 'TiAl', 'TiGa',
                       'TiGe', 'TiSb', 'TiSi', 'TiSn', 'TiZn', 'VAl', 'VGa',
                       'VGe']
        #self.compos = ['TiAl']
        # self.compos = ['NiAl', 'NiGa', 'NiGe', 'NiSi', 'NiSn', 'NiZn',
        #                'TiGa',
        #           'TiGe', 'TiSb', 'TiSi', 'TiSn', 'TiZn', 'VAl', 'VGa']
        #self.compos = ['FeAl', 'FeGa', 'FeGe', 'FeSi']
        self.compos = ['CuSi']
        self.obj = data_object

    def plt_double_MAE(self, in_dir, in_dir2, out_dir):
        """ 2014/01/28 MAEをプロット for Heusler"""
        for compo in self.compos:
            print compo
            dir_List = Bash.find_files(os.path.join(in_dir, compo),
                                                   'fixed_*')
            dir_list_total = [os.path.join(in_dir, compo, x) for x in dir_List]
            #regular
            DataA = self.obj(dir_list_total, 'POSCAR',
                             'OUTCAR_polarized', 'regular')
            for i in range(1, len(compo)):
                if compo[i].isupper():
                    ele_list = [compo[0:i], compo[i:]]

            DataA.get_enthalpy(*ele_list)
            DataA.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC001',
                          'OSZICAR_SOC100', 'MAE[100-001]')
            DataA.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC001',
                          'OSZICAR_SOC110', 'MAE[110-001]')
            DataA.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC100',
                          'OSZICAR_SOC110', 'MAE[110-100]')

            #inverse
            dir_list_total = [os.path.join(in_dir2, compo, x)
                              for x in dir_List]
            DataB = self.obj(dir_list_total, 'POSCAR', 'OUTCAR_polarized',
                             'inverse')
            DataB.get_enthalpy(*ele_list)
            DataB.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC001',
                          'OSZICAR_SOC100', 'MAE[100-001]')
            DataB.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC001',
                          'OSZICAR_SOC110', 'MAE[110-001]')
            DataB.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC100',
                          'OSZICAR_SOC110', 'MAE[110-100]')

            plt = PlotResults('Fe$_2$%s' % compo, 3)
            plt.EneMomMae_CovA(DataA.data_array, 'blue', DataA.label)
            plt.EneMomMae_CovA(DataB.data_array, 'magenta', DataB.label)
            pylab.savefig('%s/%s.eps' % (out_dir, compo))
            #pylab.show()


class PlotMAEantiP(PlotMAE):
    def __init__(self, dir_list, data_object):
        self.dir_list = ['Fe3AlB', 'Fe3GaC', 'Fe3NbN',
                         'Fe3SiB', 'Fe3TiN', 'Fe3ZrC',
                         'Fe3GaN', 'Fe3NiN', 'Fe3SiC',
                         'Fe3VB', ' Fe3ZrN', 'Fe3AlC',
                         'Fe3GeN', 'Fe3SiN', 'Fe3VC',
                         'Fe4N', 'Fe3AlN', 'Fe3PdN',
                         'Fe3SnB', 'Fe3VN', 'Fe3CoN',
                         'Fe3MnN', 'Fe3ZnC', 'Fe3CuN',
                         'Fe3NbB', 'Fe3PtC', 'Fe3TiB',
                         'Fe3ZnN', 'Fe3GaB', 'Fe3NbC',
                         'Fe3PtN', 'Fe3TiC', 'Fe3ZrB']
        self.dir_list = dir_list
        self.obj = data_object

    def plot_mae(self, in_dir, out_dir):
        for dirc in self.dir_list:
            print(dirc)
            dir_List = Bash.find_files(os.path.join(in_dir, dirc),
                                                   'fixed_*')
            dir_list_total = [os.path.join(in_dir, dirc, x) for x in dir_List]
            #regular
            DataA = self.obj(dir_list_total, 'POSCAR',
                             'OUTCAR_polarized', 'regular')
            for i in range(1, len(dirc)):
                if dirc[i].isupper():
                    ele_list = [dirc[0:i], dirc[i:]]

            DataA.get_enthalpy(*ele_list)
            DataA.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC001',
                          'OSZICAR_SOC100', 'MAE[100-001]')
            DataA.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC001',
                          'OSZICAR_SOC110', 'MAE[110-001]')
            DataA.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC100',
                          'OSZICAR_SOC110', 'MAE[110-100]')

            #inverse
            dir_list_total = [os.path.join(in_dir2, dirc, x) for x in dir_List]
            DataB = self.obj(dir_list_total, 'POSCAR', 'OUTCAR_polarized',
                             'inverse')
            DataB.get_enthalpy(*ele_list)
            DataB.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC001',
                          'OSZICAR_SOC100', 'MAE[100-001]')
            DataB.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC001',
                          'OSZICAR_SOC110', 'MAE[110-001]')
            DataB.add_mae(dir_list_total, 'POSCAR', 'OSZICAR_SOC100',
                          'OSZICAR_SOC110', 'MAE[110-100]')

            plt = PlotResults('Fe$_2$%s' % dirc, 3)
            plt.EneMomMae_CovA(DataA.data_array, 'blue', DataA.label)
            plt.EneMomMae_CovA(DataB.data_array, 'magenta', DataB.label)
            pylab.savefig('%s/%s.eps' % (out_dir, dirc))
            #pylab.show()

class AkaiKKR:
    @staticmethod
    def get_data(dir_list, output_file):
        dbox = []
        for dir in dir_list:
            path = os.path.join(dir, output_file)
            ene_line = commands.getoutput('grep \'total energy\' %s' % path)
            atom_line = commands.getoutput('grep \'natm=\' %s' % path)
            if ene_line == '':
                print '\"no output total energy\"'
                os.system('pwd')
            else:
                atom_num = int(atom_line.split()[3])
                energy = float(ene_line.split()[2]) / atom_num * 13.6058
                latt_line = commands.getoutput('grep \'brvtyp=\' %s' % path)
                a = float(latt_line.split()[2])
                cova = float(latt_line.split()[4])
                volume = ((a**3) * cova * 0.529177**3) / atom_num
                spn_line = commands.getoutput('grep \'chr,spn\' %s' % path)
                i = spn_line.split().index('chr,spn')
                mag = float(spn_line.split()[i+2]) / atom_num
                mag = math.fabs(mag)
                dbox.append([cova, energy, volume, mag])
        dbox.sort(key=lambda x: x[2])
        dtype = [('c/a', 'f8'), ('energy', 'f8'),
                 ('volume', 'f8'), ('mag', 'f8')]
        return numpy.array(map(tuple, dbox), dtype=dtype)


def trim(self, array, label, value):
        trim_list = []
        for var in array[label]:
            if var == value:
                trim_list.append(True)
            else:
                trim_list.append(False)
        trim_array = numpy.array(trim_list)
        return array[trim_array]


class VaspCovA(object):
    """
    dirctoryのlistを引数として、
    list中の全てのposcar, oszicarファイルを読んで、
    c/a, volume, energy, magのself.data_arrayを作成する
    volume, energy, magの3つは単位原子あたりに換算
    """
    def __init__(self, dir_list, poscar='POSCAR', out_file='OSZICAR'):
        self.data_array = self.get_data(dir_list, poscar, out_file)

    def get_data(self, dir_list, poscar='POSCAR', out_file='OSZICAR'):
        """
        POSCARから読み取ったc/aとvolume,
        OSZICAR/OUTCARから読み取ったenergyとmagを
        np.arrayでreturnする
        """
        data_box = []
        for dirc in dir_list:
            path_oszi = os.path.join(dirc, out_file)
            path_posc = os.path.join(dirc, poscar)
            result = self.get_data_single(path_posc, path_oszi)
            data_box.append(result)
        data_box.sort(key=lambda x: x[0])
        dtype = [('c/a', 'f8'), ('energy', 'f8'),
                 ('volume', 'f8'), ('mag', 'f8')]
        return numpy.array([tuple(x) for x in data_box], dtype=dtype)

    def get_data_single(self, poscar='POSCAR', oszicar='OSZICAR'):
        """
        src_dir中のデータを修得
        POSCARからc/aとvolumeとformula
        OSZICAR/OUTCARから読み取ったenergyとmagをreturn
        volume, energy, magの3つは単位原子あたりに換算
        """
        cova, volume, num_atoms = self.get_lattices(poscar)
        energy, mag = self.get_ene_mag(oszicar)
        result = ([cova, energy / num_atoms,
                   volume / num_atoms, mag / num_atoms])
        return result

    @staticmethod
    def get_lattices(poscar='POSCAR'):
        """
        Read cova, volume, num_atoms from POSCAR file.
        """
        posc = vaspy.Poscar(poscar)
        a, _, c = posc.get_lattice_length()
        cova = c/a
        volume = posc.get_cell_volume()
        num_atoms = sum(posc.num_atoms)
        return cova, volume, num_atoms

    @staticmethod
    def get_ene_mag(oszicar='OSZICAR'):
        """
        Read energy and magnetic momentum from OSZICAR/OUTCAR.
        """
        try:
            lines = Cabinet.read_file(oszicar)
        except IOError:
            lines = ''

        if lines[0].split()[0] == 'N':  # OSZICAR
            output = vaspy.Oszicar(oszicar)
        elif lines[0].split()[0][0:4] == 'vasp':  # OUTCAR
            output = vaspy.Outcar(oszicar)
        else:
            print("I cannot read {0}".format(oszicar))
            return None, None
        energy = output.results[-1]['energy']
        mag = output.results[-1]['mag']
        return energy, mag

    def add_mae(self, dir_list, poscar, osz1, osz2, label):
        """
        MAEの計算値をself.data_arrayに追加
        """
        array_soc1 = self.get_data(dir_list, poscar, osz1)
        array_soc2 = self.get_data(dir_list, poscar, osz2)
        mae = (array_soc2['energy'] - array_soc1['energy']) * (10 ** 6)
        dtype = (label, 'f8')
        self.data_array = self.add_data(self.data_array, mae, dtype)

    @staticmethod
    def add_data(src_array, add_row, dtype):
        """
        ラベル付きarrayに新しい一次元dataを加える⇒commopy行き??
        """
        new_dtype = src_array.dtype.descr + [dtype]
        new_array = numpy.array(map(list, src_array))
        add = []
        for i in add_row:
            add.append([i])
        new_array = numpy.append(new_array, add, 1)
        new_array = numpy.array(map(tuple, new_array), new_dtype)
        return new_array

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

    def print_array(self, *print_indexes):
        """Print data_array based on print_indexes"""
        out_lines = "\t".join(print_indexes)
        out_lines += "\n"
        for i in range(0, len(self.data_array[print_indexes[0]])):
            lines = []
            for index in print_indexes:
                lines.append(str(self.data_array[index][i]))
            out_lines += "\t".join(lines)
            out_lines += "\n"
        print(out_lines)

    def get_enthalpy(self, poscar, potcar):
        """Convert energy (eV/atom) to enthalpy (kJ/mol)"""
        elements, num_atoms = self.get_composition(poscar, potcar)
        ref_energy = 0
        total_num = sum(num_atoms)
        for element, num in itertools.izip(elements, num_atoms):
            ref_energy += Reference[element]*num
        ref_energy = ref_energy / total_num
        self.data_array['energy'] = ((self.data_array['energy'] -
                                      ref_energy) * 96.485344520851)


class HeuslerFe2YZ(VaspCovA):
    """ホイスラー用のobject"""
    def __init__(self, dir_list, poscar, oszicar, label):
        """軸比は1/√2倍する"""
        VaspCovA.__init__(self, dir_list, poscar, oszicar)
        self.data_array['c/a'] = self.data_array['c/a']/(2**0.5)
        self.label = label

    def get_enthalpy(self, compY, compZ):
        """
        POSCARから読めないのでcompY, compZを指定して計算
        OUTCARから読めるようにしたい...
        """
        ref_energy = (2*Reference['Fe'] + Reference[compY] +
                      Reference[compZ])/4
        enthalpy = ((self.data_array['energy'] -
                                      ref_energy) * 96.485344520851)
        dtype = ('enthalpy', 'f8')
        self.data_array = self.add_data(self.data_array, enthalpy, dtype)



class HeuslerX2FeZ(VaspCovA):
    """ホイスラー用X2FeZのobject"""
    def __init__(self, dir_list, poscar, oszicar, label):
        """軸比は1/√2倍する"""
        VaspCovA.__init__(self, dir_list, poscar, oszicar)
        self.data_array['c/a'] = self.data_array['c/a']/(2**0.5)
        self.label = label

    def get_enthalpy(self, compY, compZ):
        """
        POSCARから読めないのでcompY, compZを指定して計算
        OUTCARから読めるようにしたい...
        """
        ref_energy = (Reference['Fe'] + 2 * Reference[compY] +
                      Reference[compZ]) / 4
        self.data_array['energy'] = ((self.data_array['energy'] -
                                      ref_energy) * 96.485344520851)


class FittingAnalysis(object):
    """
    fitting用のobject
    fitting_analysis.pyに移行
    """
    def Stineman_interp_fit(self, x, y):
        yp = None
        xi = pylab.linspace(x[0], x[-1], 100)
        yi = pylab.stineman_interp(xi, x, y, yp)
        return xi, yi

    def localminimum(self, x, y):  # to find local minimum
        min = numpy.r_[True, y[1:] < y[:-1]] & numpy.r_[y[:-1] < y[1:], True]
        x_min = x[min]
        y_min = y[min]
        return x_min, y_min

    def Murnaghan_fit(self, volume_array, energy_array):
        v = volume_array
        e = energy_array
        vfit = numpy.linspace(min(v), max(v), 100)
        # fit a parabola to the data, y = ax^2 + bx + c
        a, b, c = pylab.polyfit(v, e, 2)

        #now here are our initial guesses.
        v0 = -b/(2*a)
        e0 = a*v0**2 + b*v0 + c
        b0 = 2*a*v0
        bP = 4

        x0 = [e0, b0, bP, v0]
        murnumpyars, ier = scipy.optimize.leastsq(self.Murnaghan_err,
                                                  x0, args=(e, v))
        return vfit, self.Murnaghan_func(murnumpyars, vfit), murnumpyars

    def Murnaghan_func(self, parameters, vol):
        """
        given a vector of parameters and volumes, return a vector of energies.
        equation From PRB 28,5480 (1983)
        """
        E0 = parameters[0]
        B0 = parameters[1]
        BP = parameters[2]
        V0 = parameters[3]
        E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)
        return E

    def Murnaghan_err(self, pars, y, x):
        #we will minimize this function
        err = y - self.Murnaghan_func(pars, x)
        return err


class PlotResults(FittingAnalysis):
    def __init__(self, title, figs_num):
        figsize_dict = [{'figsize': (600/80, 800/80/3)},
                        {'figsize': (600/80, 800/80*2/3)},
                        {'figsize': (600/80, 800/80)}]
        self.fig = pylab.figure(**figsize_dict[figs_num - 1])
        self.color_palette = solid.COLOR_PALETTE
        self.axis_dict = {'energy': 'Energy (eV/atom)',
                          'mag': 'Mag. ($\mu_B$/atom)'}
        default_font = {'font.family': 'Times New Roman'}
        pylab.rcParams.update(default_font)

        pylab.suptitle('\n %s' % title, size=26)

    def Ene_Vol_single(self, data_array, color_id, label):
        pltstyle_dict = self.mkstyle(color_id)
        pltstyle_dict.update({'linestyle': 'None'})
        data_label = label
        x, labelx = data_array['volume'], 'Volume ($\AA^3$/atom)'
        y, labely = data_array['energy'], 'Energy (eV/atom)'
        if 'self.ax1' not in dir(self):
            self.ax1 = pylab.subplot(111)
        self.ax1.plot(x, y, **pltstyle_dict)
        pylab.xlabel(labelx, fontsize=20)
        pylab.ylabel(labely, fontsize=20)
        xi, yi, Murnaghan_param = self.Murnaghan_fit(x, y)
        style = {}
        style.update({'color': pltstyle_dict['color']})
        self.add_plt(self.ax1, xi, yi, style)
        # pylab.text(0.4,0.5,"Minimum volume = %1.2f $\AA^3$"
        #            % Murnaghan_param[3], transform = self.ax1.transAxes)
        pylab.text(0.4, 0.4, 'Bulk modulus = %1.2f eV/$\AA^3$ = %1.2f GPa'
                   % (Murnaghan_param[1], Murnaghan_param[1]*160.21773),
                   transform=self.ax1.transAxes)
        return Murnaghan_param

    def mkstyle(self, color_id):
        id = color_id
        cl = self.color_palette
        style = {'color': cl[id][0], 'markeredgecolor': cl[id][0],
                 'markerfacecolor': cl[id][1]}
        style.update({'markersize': 8, 'linewidth': 1.5, 'markeredgewidth': 2,
                      'marker': 'o', 'linestyle': '-'})
        return style

    def default_frame31shareX(self, data_label, x_tpl, y1_tpl,
                              y2_tpl, y3_tpl, pltstyle_dict):
        x, labelx = x_tpl
        y1, label1 = y1_tpl
        y2, label2 = y2_tpl
        y3, label3 = y3_tpl

        pylab.subplots_adjust(hspace=0.001)
        self.ax1 = pylab.subplot(311)
        self.ax1.plot(x, y1, label=data_label, **pltstyle_dict)
        pylab.ylabel(label1, fontsize=20, horizontalalignment='left')
        pylab.yticks(fontsize=13)
        #self.ax1.yaxis.set_label_coords(-.12, 0.5)
        self.ax1.yaxis.set_label_coords(-.10, 0.25)

        pylab.ylim(-30, -10)
        self.ax1.legend()

        self.ax2 = pylab.subplot(312, sharex=self.ax1)
        self.ax2.plot(x, y2, **pltstyle_dict)
        self.ax2.yaxis.tick_right()
        pylab.ylim(0, 3)
        self.ax2.yaxis.set_label_coords(1.12, 0.5)
        pylab.ylabel(label2, fontsize=20)

        self.ax3 = pylab.subplot(313, sharex=self.ax1)
        self.ax3.plot(x, y3, **pltstyle_dict)

        xticklabels = self.ax1.get_xticklabels()+self.ax2.get_xticklabels()
        pylab.setp(xticklabels, visible=False)
        pylab.ylabel(label3, fontsize=20)
        pylab.xlabel(labelx, fontsize=20)
        pylab.xticks(fontsize=13)
        #pylab.xlim(0.6, 1.4)

    def default_frame21(self, data_label, x_tpl, y1_tpl,
                        y2_tpl, pltstyle_dict):
        x, labelx = x_tpl
        y1, label1 = y1_tpl
        y2, label2 = y2_tpl

        self.ax1 = pylab.subplot(211)
        self.ax1.plot(x, y1, label=data_label, **pltstyle_dict)
        pylab.ylabel(label1, fontsize=20, horizontalalignment='left')
        pylab.yticks(fontsize=13)
        self.ax1.yaxis.set_label_coords(-.12, 0.5)
        self.ax1.legend()

        self.ax2 = pylab.subplot(212)
        self.ax2.plot(x, y2, **pltstyle_dict)
        #self.ax2.yaxis.tick_right()
        #pylab.ylim(0,3)
        #self.ax2.yaxis.set_label_coords(1.12, 0.5)
        pylab.ylabel(label2, fontsize=20)

        #xticklabels = self.ax1.get_xticklabels()+self.ax2.get_xticklabels()
        #pylab.setp(xticklabels, visible=False)
        #pylab.ylabel(label3, fontsize = 20)
        pylab.xlabel(labelx, fontsize=20)
        pylab.xticks(fontsize=13)

    def add_plt(self, axis, x, y, style):
        axis.plot(x, y, **style)

    def alt_plt(self, axis_num, x, y, style):
        pylab.delaxes(self.ax2)
        ax = pylab.subplot(axis_num)
        ax.plot(x, y, **style)

    def EneMomVol_CovA(self, data_array, color_id, label):
        pltstyle_dict = self.mkstyle(color_id)
        data_label = label
        x = data_array['c/a'], 'c/a'
        y1 = data_array['energy'], 'Enthalpy (kJ/mol)'
        y2 = data_array['mag'], 'Mag. ($\mu_B$/atom)'
        y3 = data_array['volume'], 'Volume ($\AA^3$)'
        self.default_frame31shareX(data_label, x, y1, y2, y3, pltstyle_dict)
        #xi, yi = self.Stineman_interp_fit(x[0], y1[0])
        #x_min, y_min = self.localminimum(xi, yi)
        #style = {'marker': 'v', 'ls': 'None',
        #         'c': '#FFFF00', 'ms': 10, 'linewidth': 2.0}
        #self.add_plt(self.ax1, x_min, y_min, style)

    def Ene_Vol(self, data_array, color_id, label):
        pltstyle_dict = self.mkstyle(color_id)
        pltstyle_dict.update({'linestyle': 'None'})
        data_label = label
        x = data_array['volume'], 'Volume ($\AA^3$)'
        y1 = data_array['energy'], 'Energy (eV/atom)'
        y2 = data_array['mag'], 'Mag. ($\mu_B$/atom)'
        y3 = data_array['volume'], 'Volume ($\AA^3$)'
        self.default_frame31shareX(data_label, x, y1, y2, y3, pltstyle_dict)
        xi, yi, Murnaghan_param = self.Murnaghan_fit(x[0], y1[0])
        style = {}
        style.update({'color': pltstyle_dict['color']})
        self.add_plt(self.ax1, xi, yi, style)
        return Murnaghan_param

    def Ene_CovA(self, data_array, color_id, label):
        pltstyle_dict = self.mkstyle(color_id)
        pltstyle_dict.update({'linestyle': 'None'})
        data_label = label
        x = data_array['volume'], 'c/a'
        y1 = data_array['energy'], 'Energy (eV/atom)'
        y2 = data_array['mag'], 'Mag. ($\mu_B$/atom)'
        self.default_frame21(data_label, x, y1, y2, pltstyle_dict)
        xi, yi, Murnaghan_param = self.Murnaghan_fit(x[0], y1[0])
        style = {}
        style.update({'color': pltstyle_dict['color']})
        self.add_plt(self.ax1, xi, yi, style)
        return Murnaghan_param

    def EneMomMae_CovA(self, data_array, color_id, label):
        pltstyle_dict = self.mkstyle(color_id)
        data_label = label
        x = data_array['c/a'], 'c/a'
        y1 = data_array['energy'], 'Energy (kJ/mol)'
        y2 = data_array['mag'], 'Mag. ($\mu_B$/atom)'
        y3 = data_array['MAE[100-001]'], 'MAE ($\mu$eV/atom)'
        #y3_2 = data_array['MAE[110-001]'], ''
        self.default_frame31shareX(data_label, x, y1, y2, y3, pltstyle_dict)
        xi, yi = self.Stineman_interp_fit(x[0], y1[0])
        x_min, y_min = self.localminimum(xi, yi)
        style = {'marker': 'o', 'ls': 'None', 'c': 'blue',
                 'ms': 8, 'linewidth': 2.0}
        #self.add_plt(self.ax3, x[0], y3_2[0], style)
        style = {'marker': 'v', 'ls': 'None', 'c': '#FFFF00',
                 'ms': 10, 'linewidth': 2.0}
        #self.add_plt(self.ax1, x_min, y_min, style)

    def chLabelY(self, num, new_label):
        num = 310 + num
        pylab.subplot(num)
        pylab.ylabel(new_label)

    def choiceLower(self, row1, label1, row2, label2):
        new_row = []
        for i in range(0, len(row1)):
            if row1[i] < row2[i]:
                new_row.append(row1[i])
            else:
                new_row.append(row2[i])
        return new_row

def readFile(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close
    return lines


def writeFile(fname, lout):
    out = open(fname, 'w')
    for line in lout:
        out.write(line)
    out.close

if __name__ == "__main__":
    main()
