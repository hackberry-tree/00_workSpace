#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
sprkkr の rws を修正する
空間群を考慮していないのでユニットセル中のサイト数が実際と違う場合がある
その場合は corr_site を使って補正する
補正が必要ない場合は --site は使わない
"""
import argparse
import glob
from sprkkr import Pot


def main():
    """main"""
    parser = argparse.ArgumentParser(description='correct rws of pot')
    parser.add_argument('--site', dest='site', type=int, nargs=1,
                        action='store', default=[1])
    args = parser.parse_args()
    correct_rws(args.site[0])


def correct_rws(corr_site):
    """
    rmt, rwsを変更する
    """
    fname = glob.glob('*.pot')
    print(fname)
    pot = Pot.from_file(fname[0])
    pot.correct_rws(corr_site)
    wfile = fname[0] + "_rws"
    pot.write_file(wfile)


if __name__ == '__main__':
    main()
