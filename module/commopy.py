#!/opt/local/bin/python
# -*- coding: utf-8 -*-
"""
My trivial python module
"""
import datetime
import os
import socket
import copy
import glob
import shutil
import errno
import math
import subprocess
import numpy as np


class LinkedList(object):
    """連結リスト構造"""

    class Cell(object): #pylint: disable=R0903
        """
        セル用のクラス 他では使用ししない
        data: data
        link: 次の参照セルへのlink
        """
        def __init__(self, data, link=None):
            self.data = data
            self.link = link

    def __init__(self, *items):
        self.top = LinkedList.Cell(None)   # ヘッダセル
        for item in reversed(items):
            self.insert(0, item)

    def _nth(self, n): #pylint: disable=C0103
        """n番目のセルを求める"""
        i = -1
        cell = self.top
        while cell is not None:
            if i == n:
                return cell
            i += 1
            cell = cell.link
        return None

    def at(self, n):
        """n番目の要素を求める"""
        cell = self._nth(n)
        if cell is not None:
            return cell.data
        return None

    def insert(self, n, data):
        """データの挿入"""
        cell = self._nth(n - 1)
        if cell is not None:
            cell.link = LinkedList.Cell(data, cell.link)
            return data
        return None

    def delete(self, n):
        """データの削除"""
        cell = self._nth(n - 1)
        if cell is not None and cell.link is not None:
            data = cell.link.data
            cell.link = cell.link.link
            return data
        return None

class DataBox(object):
    """
    labeled list data
    list-typeのdata格納class
    dataの呼び出し時にarrayに変換して演算処理が可能
    """
    def __init__(self, data_list):
        """リストデータを作成"""
        self.data = data_list
        self.output_keys = []

    def update(self, key, add_data):
        """
        self.dataにデータを追加
        keyが無い場合新しく作成
        """
        if len(self.data) != len(add_data):
            print("Dimension of data is different between self.data ",
                  "and add_data !!!")
        for value, src_data in zip(add_data, self.data):
            src_data.update({key: value})

    def set_order(self):
        """
        順番のラベルを入れる
        横軸がtextのgraphを出力する際に利用
        """
        order = [i for i in range(0, len(self.data))]
        self.update('order', order)

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

    def trim_data(self, distinct, value):
        """
        distinct(list or array)がvalueと等しいdataだけを抜き出す
        """
        if len(distinct) != len(self.data):
            print("Size of data is differnt !!")
        return [self.data[i] for i in range(0, len(self.data))
                if distinct[i] == value]

    def separate_data(self, key):
        """
        self.dataをself.data[key]が同じ値を持つデータ毎に分割
        DataBox形式のデータをもつlistをreturn
        """
        values = sorted(set(self[key]))
        separated_data = []
        for val in values:
            data = [x for x in self.data if x[key] == val]
            separated_data.append(DataBox(data))
        return separated_data

    def __str__(self):
        """
        self.output_keysで指定した項目を出力
        指定無ければ、全てを出力
        """
        keys = self.output_keys
        if not keys:
            keys = self.data[0].keys()  # output_keysが指定なければ全てを出力

        # dataがdictの場合ラベルを細分化
        types_dict = [x for x in keys if type(self.data[0][x]) is dict]
        out_lines = "\t".join(keys) + "\n"
        for key in types_dict:
            head = "".join(key.split('_dict')) + "_"
            label = [head + str(x) for x in self.data[0][key]]
            alt_label = "\t".join(label)
            out_lines = out_lines.replace(key, alt_label)
        for data in self.data:
            line = []
            for key in keys:
                if type(data[key]) is dict:
                    data[key] = "\t".join([str(x) for x in data[key].values()])
                line.append(str(data[key]))
            out_lines += "\t".join(line) + "\n"
        return out_lines

    def to_array(self):
        """
        arrayに変換
        """
        #pylint: disable=E1101
        return np.array([[x[key] for key in self.output_keys] \
                        for x in self.data])

    def to_list(self):
        """
        output_keysのみを要素にもつlistをreturn
        """
        out = []
        for data in self.data:
            single = {key: data[key] for key in self.output_keys}
            out.append(single)
        return out

# class Convert(object):
#     """
#     何か変換する
#     """

#     @staticmethod
#     def strings(attributes):
#         """
#         変数名を文字列のlistに変える
#         多重参照の場合 任意性(len(list) > 1)があるので注意
#         errorを特定するのに使える程度
#         """
#         name = [key for key, value in globals().items()
#                              if id(attributes) == id(value)]
#         return name


class Cabinet(object):
    """
    Handling files

    attributes: read_file(fname), write_file(fname, lout), reserve_file(file)
    """

    @staticmethod
    def read_file(fname):
        """
        read_file(fname) read 'fname' file and return list lines.
        """
        with open(fname, 'r') as fread:
            read_lines = fread.readlines()
        return read_lines

    @staticmethod
    def read_file_1line(fname):
        """
        read_file(fname) read 'fname' file and return list lines.
        """
        with open(fname, 'r') as fread:
            read = fread.read()
        return read

    @staticmethod
    def write_file(fname, lines):
        """
        write_file(fname, lines) write 'lines' into 'fname' file.
        """
        with open(fname, 'w') as fwrite:
            for line in lines:
                fwrite.write(line)

    @staticmethod
    def append_file(fname, lines):
        """
        Append lines into fname file.
        """
        with open(fname, 'a') as fappend:
            for line in lines:
                fappend.write(line)

    @staticmethod
    def reserve_file(fname):
        """
        If old file exit, it is reserved with adding time stamp in file name.
        """
        if glob.glob(fname):
            date = datetime.datetime.today()
            stamp = date.strftime("_%y%b%d_%H%M%S")
            alt_name = fname + stamp
            Bash.move(fname, alt_name)
            print("Old '{0}' file is saved as '{1}'".format(fname, alt_name))

    @staticmethod
    def conv_line2lines(line):
        """
        改行でsplitしてそれぞれの末尾に'\n'を追加したlistをreturnする
        read_fileの出力と対応したものを得る
        主にstr(obj)の出力を変換したい時に使用
        """
        lines = [x + '\n' for x in line.split('\n')]
        return lines

    @staticmethod
    def conv_lines2array(lines_list, dtype=float, split=None):
        """
        read_fileから読み込んだスペース区切りか
        タブ区切りのlinesからなるlistをarrayに変換してreturnする
        dtypeを変更することで出力データのtypeを指定
        区切りの種類はsplitを指定することで変更可
        """
        box = []
        for line in lines_list:
            box.append([dtype(x) for x in line.split(split)])
        return np.array(box) #pylint: disable=E1101

    @staticmethod
    def conv_str(string):
        """
        strを変換
        int > float (> str)の順に変換を試してreturn
        """
        try:
            return int(string)
        except ValueError:
            pass
        try:
            return float(string)
        except ValueError:
            return string

    @staticmethod
    def make_list_run(path_list, run_file):
        """
        next.py用のlist_run.txtを作成
        path_listとrun_fileには相対pathを入力
        outputは絶対パスで作成
        """
        lines = ""
        for path in path_list:
            path_dir = os.path.abspath(path)
            run_file = os.path.abspath(run_file)
            lines += "{0}    {1}\n".format(path_dir, run_file)
        Cabinet.write_file('list_run.txt', lines)


class TrialRun(object):
    """
    When the condition is satisfied,
    perform the registered functions.
    """
    HOST = socket.gethostname()
    INITIAL_LINES = ("\n%======="
                     "This is TEST_RUN @ {0}"
                     "=======%\n\n".format(HOST))
    END_LINES = ("\n%====="
                 "TEST_RUN is finished @ {0}"
                 "======%\n".format(HOST))

    @classmethod
    def execute(cls, func):
        """
        Defualt execute function
        """
        print(cls.INITIAL_LINES)
        func()
        print(cls.END_LINES)
        exit()

    @classmethod
    def is_workspace(cls, dirc, func):
        """
        When the current directory is workspace,
        perform a registered function (test run).
        """
        if os.path.abspath(dirc) == ("/Users/enoki/Documents/01_ResearchData/"
                                     "Calculations/99_python/00_workSpace"):
            cls.execute(func)


class Bash(object):
    """bashでのコマンド関連"""
    @staticmethod
    def execute(cmd):
        """cmdをbash上で実行"""
        try:
            retcode = subprocess.call(cmd, shell=True)
            if retcode < 0:
                print("Stopped \'{0}\' process from signal"
                      " {1}".format(cmd, -retcode))
            else:
                print("\'{0}\' process is normally finished"
                      " {1}".format(cmd, retcode))
        except OSError as err:
            print("Faile to execute \'{0}\' process:"
                  " {1}".format(cmd, err))

    @staticmethod
    def back_quote(*cmd):
        """cmdの実行結果をreturn"""
        output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
        return output

    @staticmethod
    def find_files(fname):
        """
        Find file in path
        Enable to use wild card "*"
        globのまま使うので十分か...
        一応 basenameと切り分けてファイル名のみを出力する
        """
        dir_list = [os.path.basename(x) for x in glob.glob(fname)]
        return dir_list

    @staticmethod
    def copy_dir(src_dir, dst_dir):
        """
        Copy directory to dst_dir.
        If dst_dir exit. Do Nothing.
        重複以外残して上書きする記述もありかもしれない、保留
        ファイルのコピーに関しては既存のファイルがあっても上書きするので、
        shutil.copyfile(src, dst)をそのまま使う
        """
        if os.path.exists(dst_dir):
            print("\"{0}\" already exist.\n"
                  "Do Nothing...\n".format(dst_dir))
            return
        shutil.copytree(src_dir, dst_dir)

    @staticmethod
    def mkdir(path):
        """
        Alternative command of 'mkdir -p'
        At Python3.3, we enable to use os.makedirs(path, exist_ok=True)
        """
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    @staticmethod
    def move(src_path, dst_path):
        """
        'mv' command.
        If old file exists, old file will be over_writed using os.rename.
        dst_pathにファイルがあって、src_pathが無い場合、
        ファイルが消えてしまうようなので例外処理で防止する
        """
        try:
            shutil.move(src_path, dst_path)
        except shutil.Error as err:
            if not os.path.exists(src_path):
                print("Error: {0} dose not exists !".format(src_path))
                exit()
            dst_file = str(err).split(' \'')[-1].split('\' ')[0]
            print(dst_file)
            os.rename(src_path, dst_file)


class Compare(object):
    """比較する"""
    @staticmethod
    def dict(dict1, dict2):
        """
        Compare two dictionary.
        Out put : diffrence between the two dictionary.
        """
        key_list = list(dict1.keys()) + list(dict2.keys())
        comp1 = {}
        for key in set(key_list):
            comp1.update({key: None})
        comp2 = copy.deepcopy(comp1)
        comp1.update(dict1)
        comp2.update(dict2)
        count = 0.
        correct = 0.
        for key in set(key_list):
            if comp1[key] == comp2[key]:
                correct += 1
            else:
                print("{0} is {1} and {2}".format(key, comp1[key], comp2[key]))
            count += 1
        frac = correct / count * 100
        print("Two dict have {0:.1f}% same contents".format(frac))
        return frac


class Array(object):
    """numpy array関連"""
    @staticmethod
    def add_data(src_array, add_row, dtype):
        """
        ラベル付きarrayに新しい一次元dataを加える⇒commopy行き??
        """
        new_dtype = src_array.dtype.descr + [dtype]
        new_array = np.array([list(x) for x in src_array]) #pylint: disable=E1101
        add = []
        for i in add_row:
            add.append([i])
        new_array = np.append(new_array, add, 1) #pylint: disable=E1101
        new_array = np.array([tuple(x) for x in new_array], dtype=new_dtype) #pylint: disable=E1101
        return new_array

    @staticmethod
    def extract(src_array, *labels):
        """
        labelsのデータのみを抽出してarrayをreturn
        """
        new_dtype = []
        for key, dtype in src_array.dtype.descr:
            if key in labels:
                new_dtype.append((key, dtype))
        new_data = [tuple(src_array[key][x] for key in labels)
                    for x in range(src_array.shape[0])]
        return np.array(new_data, dtype=new_dtype) #pylint: disable=E1101

    @staticmethod
    def trim(array, label, value):
        """
        labelsのデータがvalueと等しいデータのみを抽出して
        arrayをreturnする
        数値のみに対応
        """
        trim_list = []
        for var in array[label]:
            if (var - value) ** 2 < 1e-16:
                trim_list.append(True)
            else:
                trim_list.append(False)
        trim_array = np.array(trim_list) #pylint: disable=E1101
        return array[trim_array]

    @staticmethod
    def trim_bool(array, label, value):
        """
        labelsのデータがvalueと等しいデータのみを抽出して
        boolをreturnする
        数値のみに対応
        """
        trim_list = []
        for var in array[label]:
            if (var - value) ** 2 < 1e-16:
                trim_list.append(True)
            else:
                trim_list.append(False)
        return np.array(trim_list) #pylint: disable=E1101


    @staticmethod
    def frange_np(start, stop, num_points):
        """
        rangeの小数版
        ポイントの数を指定
        """
        out = []
        stp = (stop - start) / (np - 1)
        for i in range(num_points):
            out.append(start + stp * i)
        return out

    @staticmethod
    def frange_stp(start, stop, step):
        """
        rangeの小数版
        step間隔を指定する
        """
        out = []
        num_points = int((stop - start) / step + 1)
        for i in range(num_points):
            out.append(start + step * i)
        return out


class Vector(object):
    """Calclulation for vectors"""
    @staticmethod
    def get_angle(vector_a, vector_b):
        """
        Calculate angle of two vectors.
        """
        len_a, len_b = (np.linalg.norm(x) for x in (vector_a, vector_b)) #pylint: disable=E1101
        cosine = np.dot(vector_a, vector_b) / len_a / len_b #pylint: disable=E1101
        degree = math.acos(cosine) / (2 * math.pi) * 360
        return degree

    @staticmethod
    def get_volume(vector_a, vector_b, vector_c):
        """
        Calculate hexahedral volume.
        """
        matrix = [vector_a, vector_b, vector_c]
        volume = np.linalg.det(matrix) #pylint: disable=E1101
        math.fabs(volume)
        return volume
