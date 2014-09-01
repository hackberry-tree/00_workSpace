#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
list_runを作成
"""
import os
import glob
import argparse
from commopy import Cabinet

__date__ = "Sep 1 2014"

def main():
    """main"""
    execute_parser()

def execute_parser():
    """
    parser

    usage:
    fileで検索 outputのpathにはfileは含まない
    make_list_run.py '*/*/INCAR' --file --run 'calc.sh'

    dirで検索
    make_list_run.py '*/*/' --run 'calc.sh'

    複数のkeyを参照する場合 ','で区切る
    (重複のpathはsetで消される)
    make_list_run.py '*/*/INCAR, */*/POTCAR' --file --run 'calc.sh'
    """
    parser = argparse.ArgumentParser(description='make list_run.txt')
    parser.add_argument('keys', type=str, nargs='?')
    parser.add_argument('--file', dest='is_file', default=False,
                        const=True, action='store_const')
    parser.add_argument('--run', dest='run', type=str, nargs='?',
                        default='run.sh', const=str)
    args = parser.parse_args()
    keys = args.keys.split(',')
    make_list(keys, is_file=args.is_file, run_file=args.run)

def make_list(keys, is_file=False, run_file='run.sh'):
    """
    current directory中のkeyで引っ掛かるpathを
    list_run.txtに書き出す
    """
    path_list = []
    for key in keys:
        path_list += glob.glob(key)
    path_list = sorted(set(path_list))
    if is_file:
        path_list = [os.path.dirname(x) for x in path_list]
    Cabinet.make_list_run(path_list, run_file)


if __name__ == '__main__':
    main()
