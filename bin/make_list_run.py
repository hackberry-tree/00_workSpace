#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
list_run.txtを作成
ToDo:argparseを使って実行する形式にする
"""
import glob
import argparse
from commopy import Cabinet

def execute_parser():
    """
    parser
    keyは複数設定可能

    usage
    make_list_run.py '*/*/INCAR' --run 'calc.sh'
    """
    parser = argparse.ArgumentParser(description='make list_run.txt')
    parser.add_argument('keys', type=str, nargs='+')
    parser.add_argument('--run', dest='run', type=str, nargs='?',
                        default='run.sh', const=str)
    args = parser.parse_args()
    make_list(args.keys, run_file=args.run)

def main():
    """main"""
    execute_parser()


def make_list(keys, run_file='run.sh'):
    """
    current directory中のkeyで引っ掛かるpathを
    list_run.txtに書き出す
    """
    path_list = []
    for key in keys:
        path_list += glob.glob(key)
    path_list = sorted(set(path_list))
    Cabinet.make_list_run(path_list, run_file)


if __name__ == '__main__':
    main()
