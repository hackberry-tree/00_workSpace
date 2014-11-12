#!/opt/local/bin/python
# -*- coding: utf-8 -*-
"""
This module dose not work in python 3
only use python 2.7+
use pymatgen
"""
from pymatgen.matproj.rest import MPRester #pylint: disable=F0401
from pymatgen.phasediagram.pdmaker import PhaseDiagram #pylint: disable=F0401
from pymatgen.phasediagram.plotter import PDPlotter #pylint: disable=F0401
from pymatgen.io.cifio import CifParser
from pymatgen.io.vaspio import Poscar

def test():
    """for test"""
    #This initializes the REST adaptor. Put your own API key in.
    mpr = MPRester("WTxsDhRV7g2Mcbqw")
    #Entries are the basic unit for thermodynamic
    #and other analyses in pymatgen.
    #This gets all entries belonging to the Ca-C-O system.
    composition = ['Fe', 'Ni', 'Si']
    entries = mpr.get_entries_in_chemsys(composition)

    # entryの内容を確認
    entry = entries[0]
    print(entry.as_dict())
    print(entry.as_dict()['data']['volume']/6)

    #With entries, you can do many sophisticated analyses,
    #like creating phase diagrams.
    pd = PhaseDiagram(entries) #pylint: disable=C0103

    #print(dir(pd))
    #print(pd.get_form_energy_per_atom(list(pd.stable_entries)[0]))

    plotter = PDPlotter(pd)
    #plotter.show()

    #get chmical fomular
    #print(entries[0].data['cif'])  # cif data
    #cif = [CifParser.from_string(x.data['cif']) for x in entries]
    #cif_str = [x.get_structures()[0] for x in cif]
    #poscar = Poscar(cif_str[1])
    #print(poscar)

    #print(entries[0].to_dict['data'])
    #print_data(entries)
    #print_data2(entries, composition)

def print_data(entries):
    line = "formula\tenthalpy\n"
    for entry in entries:
        formula = entry.to_dict['data']['unit_cell_formula']
        formation_e = entry.to_dict['data']['formation_energy_per_atom']
        if formation_e <= 0:
            for elements, index in formula.items():
                line += elements + str(int(index))
            line += "\t" + str(formation_e) + "\n"
    print(line)

def print_data2(entries, composition):
    line = "{0}\tenthalpy\n".format("\t".join(composition))
    for entry in entries:
        formula = entry.to_dict['data']['unit_cell_formula']
        formation_e = entry.to_dict['data']['formation_energy_per_atom']
        if formation_e <= 0:
            sum_atoms = 0
            for elements, index in formula.items():
                sum_atoms += index
            for element in composition:
                try:
                    line += str(formula[element] / sum_atoms) + "\t"
                except KeyError:
                    line += "0" + "\t"
            line += "\t" + str(formation_e) + "\n"
    print(line)

if __name__ == "__main__":
    test()
