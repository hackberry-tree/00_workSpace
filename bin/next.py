#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
import os


def writeFile(fname, lout):
    out = open(fname, 'w')
    for line in lout:
        out.write(line)
    out.close


def readFile(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close
    return lines


def altLine(lines, oldword, newword):
    """
    cd $PBS_O_WORKDIRをpathに置き換える
    """
    if lines.count(oldword) != 0:
        pos = lines.index(oldword)
        lines[pos] = newword
    else:
        pass


def altRunFile(run_file, curdir, outname):
    """
    cd $PBS_O_WORKDIRをpathに置き換えて実行
    """
    oldword = 'cd $PBS_O_WORKDIR\n'
    newword = 'cd %s\n' % os.path.join(curdir, run_file[0])
    path = os.path.join(curdir, *run_file)
    lines = readFile(path)
    altLine(lines, oldword, newword)
    lines.append('cd %s\n' % curdir)
    lines.append('next.py\n')
    writeFile(outname, lines)


def main():
    lines_list = readFile('list_run.txt')
    pos = 0
    for line in lines_list:
        if '#' in line.split()[0]:
            pos = pos + 1
        else:
            next = line.split()
            break
    lines_list[pos] = '# %s' % (lines_list[pos])
    writeFile('list_run.txt', lines_list)
    curDir = os.getcwd()
    print(os.path.join(curDir, next[0]))
    altRunFile(next, curDir, '.next_run.sh')
    os.system('sh .next_run.sh')

if __name__ == "__main__":
    main()
