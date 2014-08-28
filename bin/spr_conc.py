#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""runlistの生成"""
from makeSeries import Combinatorial
from commopy import Cabinet

def main():
    """main"""
    change_conc()


def change_conc():
    """
    concを変更する
    """
    list1 = [str(i*0.1) for i in range(0, 11)]
    list2 = ['001', '100']
    list3 = ['0.950', '1.000', '1.050', '1.100', '1.150',
             '1.200', '1.250', '1.300', '1.350', '1.400',
             '1.450', '1.500']
    combi = Combinatorial(list1, list2, list3)
#    print(combi.compositions)
    comp_list = [x['elements']
                 for x in combi.compositions]
    lines = ""
    for comp in comp_list:
        lines += "/".join(comp) + " run.sh\n"
    print(lines)
    Cabinet.write_file('list_run.txt', lines)

    # path =
    # Bash.mkdir('combi')
    # elem_list = ['Ni', 'Mn', 'Cu', 'Zn', 'Cr']
    # comp_list = [['O', 'Fe', x] for x in elem_list]
    # series = series_vasp.Produce('POSCAR', 'combi')
    # series.set_elements(comp_list)
    # series.make_files()
    # series.append_list_run('run_combi.sh')

if __name__ == '__main__':
    main()
