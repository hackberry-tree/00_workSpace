#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
convex_hull を作成
Ground State を探す
"""
import copy
import math
import pylab
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from commopy import Cabinet
from mpl_toolkits.mplot3d import proj3d

default_font = {'font.family': 'Times New Roman'}
pylab.rcParams.update(default_font)

def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, a, b],
                     [0, 0, 0, zback]])


proj3d.persp_transformation = orthogonal_proj

def draw_convex_hull(ax, initial_base, not_bases, elements,
                     zlim='auto', zlabel=True, **kwargs):
    """
    convex_hullを描く
    """
    bases = FindGS.collect_base_triangles(initial_base, not_bases)
    data = np.array(initial_base + not_bases)

    pt3d = PlotTriangularCoord(ax)

    pt3d.plt_dot(data, **kwargs)
    pt3d.outline_base_plot(bases)
    if zlim == 'auto':
        zmax = math.floor(max([x[2] for x in data]) / 10) * 10 + 10
        zmin = math.ceil(min([x[2] for x in data]) / 10) * 10 - 10
        zlim = [zmin, zmax]

    pt3d.make_flame(elements, zlim)
    return ax

class PlotTriangularCoord(object):
    """
    三角座標での3Dプロット
    """
    def __init__(self, ax):
        self.ax = ax

    def plt_line(self, X, Y, Z, color='gray'):
        self.ax.plot3D(np.ravel(X), np.ravel(Y), np.ravel(Z), color=color)

    def make_flame(self, elements, zlim, zlabel=True):
        for i in range(0, 10):
            X_orig = np.array([1 - i*0.1, 0])
            Y_orig = np.array([0, 1 - i*0.1])
            X, Y = alt_triangle(X_orig, Y_orig)
            Z = np.array([0, 0])
            self.plt_line(X, Y, Z)
        for i in range(0, 10):
            X_orig = np.array([i*0.1, i*0.1])
            Y_orig = np.array([0, 1 - i*0.1])
            X, Y = alt_triangle(X_orig, Y_orig)
            Z = np.array([0, 0])
            self.plt_line(X, Y, Z)
        for i in range(0, 10):
            X_orig = np.array([0, 1 - i*0.1])
            Y_orig = np.array([i*0.1, i*0.1])
            X, Y = alt_triangle(X_orig, Y_orig)
            Z = np.array([0, 0])
            self.plt_line(X, Y, Z)
        for i in range(0, 3):
            x = [1, 0, 0]
            y = [0, 1, 0]
            X_orig = np.array([x[i], x[i]])
            Y_orig = np.array([y[i], y[i]])
            X, Y = alt_triangle(X_orig, Y_orig)
            Z = [0, zlim[0]]
            self.plt_line(X, Y, Z)
        for i in range(0, 3):
            x = [1, 0, 0, 1]
            y = [0, 1, 0, 0]
            X_orig = np.array([x[i], x[i+1]])
            Y_orig = np.array([y[i], y[i+1]])
            X, Y = alt_triangle(X_orig, Y_orig)
            Z = [zlim[0], zlim[0]]
            self.plt_line(X, Y, Z)
        x = [-0.01, 0.01]
        y = [-0.01, 0.01]
        z = [zlim[0]/2, zlim[0]/2]
        self.plt_line(x, y, z)


        self.ax.set_zlim(*zlim)
        # 回転中心が重心になるようにxlim & ylimを設定
        self.ax.set_ylim([-math.sqrt(3)/6., math.sqrt(3)/2.])
        self.ax.set_xlim([-math.sqrt(3)/3.+1/2., +math.sqrt(3)/3.+1/2.])

        self.ax.axis("off")
        #self.ax.text(0, 0, zlim[0]/2., 'Energy')
        self.ax.text(0.5, math.sqrt(3)/2+0.1, 0, elements[0], size=18)
        #center = np.array([[1/3., 1/3., 0]])
        #self.plt_dot(center)
        self.ax.text(-0.1, -0.05, 0, elements[1], size=18)
        self.ax.text(1.1, -0.05, 0, elements[2], size=18)
        #self.ax.text(0, 0, zlim[1], zlim[1], size=18)
        if zlabel:
            self.ax.text(-0.2, -0.1, zlim[0], zlim[0], size=18)
            self.ax.text(-0.2, -0.1, zlim[0]/2, int(zlim[0]/2), size=18)

    def outline_base_plot(self, bases, **kwargs):
        """
        basesを三角形状に結んでプロット
        """
        for three_points in bases:
            triangle = np.array(three_points + [three_points[0]])
            x, y = alt_triangle(triangle[:, 0], triangle[:, 1])
            self.plt_line(x, y, triangle[:, 2], **kwargs)

    def plt_dot(self, data, **kwargs):
        x, y, z = (data[:, i] for i in range(0, 3))
        tri_x, tri_y = alt_triangle(x, y)
        self.ax.scatter3D(np.ravel(tri_x), np.ravel(tri_y), np.ravel(z),
                          **kwargs)

class iCVM_energies(object):
    """
    iCVM (i-s 三元系) の energies.txt を取り扱う
    """
    VAC_SITE = 'b'
    def __init__(self, phase, num_atoms, energy, sites):
        self.num_atoms = num_atoms
        self.energy = energy
        self.sites = sites
        vac = num_atoms[:, {'d': 3, 'b': 1}[self.VAC_SITE]]
        self.total_atom = sites - vac
        self.energy_per_atom = energy / self.total_atom
        self.fract_per_atom = num_atoms / self.total_atom.reshape(-1, 1)
        self.fract_per_site = num_atoms / sites.reshape(-1, 1)
        self.enthalpy = self.get_enthalpy()

    @classmethod
    def from_file(cls, fname):
        """
        energies.txt から読み込む
        """
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        phase = [x.split('*')[0] for x in lines]
        elem_a = np.array(
            [int(x.split()[0].split('A')[1].split('B')[0]) for x in lines])
        elem_b = np.array(
            [int(x.split()[0].split('B')[1].split('C')[0]) for x in lines])
        elem_c = np.array(
            [int(x.split()[0].split('C')[1].split('D')[0]) for x in lines])
        elem_d = np.array([int(x.split()[0].split('D')[1]) for x in lines])
        energy = np.array([float(x.split()[2]) for x in lines])
        sites = np.array([int(x.split()[3]) for x in lines])
        num_atoms = np.c_[elem_a, elem_b, elem_c, elem_d]
        if not (sites == num_atoms.sum(axis=-1)).all():
            print('error')
        return cls(phase, num_atoms, energy, sites)

    def get_enthalpy(self):
        end1 = \
            (self.fract_per_atom[:, [0, 2, 3]] == [0, 1, 0]).prod(axis=-1) == 1
        end2 = \
            (self.fract_per_atom[:, [0, 2, 3]] == [0, 0, 1]).prod(axis=-1) == 1
        end3 = \
            (self.fract_per_atom[:, [0, 2, 3]] == [3/4, 1/4, 0]).prod(axis=-1) == 1
        end4 = \
            (self.fract_per_atom[:, [0, 2, 3]] == [3/4, 0, 1/4]).prod(axis=-1) == 1

        ref_i = (self.energy_per_atom[end3] + self.energy_per_atom[end4]) / 2

        ref1 = self.energy_per_atom[end1]
        ref2 = self.energy_per_atom[end2]
        ref = (ref1 * self.fract_per_atom[:, 2] +
               ref2 * self.fract_per_atom[:, 3] +
               ref_i * self.fract_per_atom[:, 0])
        return self.energy_per_atom - ref

    def get_enthalpy_per_site(self):
        end1 = \
            (self.fract_per_atom[:, [0, 2, 3]] == [0, 1, 0]).prod(axis=-1) == 1
        end2 = \
            (self.fract_per_atom[:, [0, 2, 3]] == [0, 0, 1]).prod(axis=-1) == 1
        end3 = \
            (self.fract_per_atom[:, [0, 2, 3]] == [3/4, 1/4, 0]).prod(axis=-1) == 1
        end4 = \
            (self.fract_per_atom[:, [0, 2, 3]] == [3/4, 0, 1/4]).prod(axis=-1) == 1

        ref_i = (self.energy_per_atom[end3] + self.energy_per_atom[end4]) / 2

        ref1 = self.energy_per_atom[end1]
        ref2 = self.energy_per_atom[end2]
        ref = (ref1 * self.fract_per_atom[:, 2] +
               ref2 * self.fract_per_atom[:, 3] +
               ref_i * self.fract_per_atom[:, 0])
        return self.energy_per_atom - ref


def alt_triangle(array_x, array_y):
    """
    三角図の座標に変換する
    (1, 0, 0) -> (0.5, √3/2)
    (0, 1, 0) -> (0, 0)
    (0, 0, 1) -> (1, 0)
    """
    triangle_x = (1 - 0.5 * array_x - array_y)
    triangle_y = (math.sqrt(3) * 0.5 * array_x)
    return triangle_x, triangle_y


def read_float(fname):
    """
    ファイルをfloatデータとして読み込み
    float化できない文字がある行はスキップ
    """
    stock = []
    rows = Cabinet.read_file(fname)
    for row in rows:
        row = row.split()
        try:
            row = [float(x) for x in row]  # 全てをfloat化
            stock.append(row)
        except ValueError:
            pass
    return stock
#######################################


###数値データ列の選択 入力:選択list###
def selData(data, selList):
    newData = []
    for val in data:
        element = []
        for num in selList:
            element.append(val[num])
        newData.append(element)
    return newData
#################################

class ConvertData(object):
    """
    Dataを変換する
    """
    @staticmethod
    def get_enthalpy(base, notbase):
        """enthalpyに変換"""
        ref = [[x[0]*x[2], x[1]*x[2], (1-x[0]-x[1])*x[2]] for x in base]
        ref = [ref[0][0], ref[1][1], ref[2][2]]
        base_ent = [[x[0], x[1], (x[2]-x[0]*ref[0]-x[1]*ref[1]-(1-x[0]-x[1])*ref[2])*96.485344520851] for x in base]
        notbase_ent = [[x[0], x[1], (x[2]-x[0]*ref[0]-x[1]*ref[1]-(1-x[0]-x[1])*ref[2])*96.485344520851] for x in notbase]
        return base_ent, notbase_ent


class FindGS(object):
    """
    Ground Statesを抽出する
    np.arrayの形式でデータを受け取るとエラーがでる
    とりあえずリスト化してデータを渡す
    initial_base = [list(x) for x in data[0:3, [0, 1, 3]]]
    """

    def __init__(self, data):
        pass

    @staticmethod
    def vector(pt2, pt1):
        """
        二点point1, point2に対するベクトル(array)を作成してリターン
        """
        return np.array(pt2) - np.array(pt1)

    @staticmethod
    def coeff_synthetic_vect(vec_a, vec_b, vec_c):
        """
        2つの2次元ベクトル(vec_a, vec_b)を合成してvCを表現
        t*vec_a + u*vec_b = vec_cとなるt, uを計算
        """
        norm = 1 / (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0])
        t = norm * (vec_c[0] * vec_b[1] - vec_c[1] * vec_b[0])
        u = norm * (-vec_c[0] * vec_a[1] + vec_c[1] * vec_a[0])
        return [t, u]

    @classmethod
    def relative_energy(cls, base, point):
        """
        与えられた3点(base)のenergy surfaceに対するpointの相対エネルギーを算出
        format of base: [A1, A2, A3]=[[x1, y1, e1],[x2, y2, e2],...]
        """
        base_va = cls.vector(base[1], base[0])
        base_vb = cls.vector(base[2], base[0])
        vec_point = cls.vector(point, base[0])
        coeff = cls.coeff_synthetic_vect(base_va, base_vb, vec_point)
        ref_energy = coeff[0] * base_va[2] + coeff[1] * base_vb[2]
        return vec_point[2] - ref_energy

    @staticmethod
    def is_parallel(vec1, vec2):
        """
        二つのvectorが平行か判定
        """
        if (vec1[0] * vec2[1] - vec1[1] * vec2[0]) ** 2 < 1e-16:
            return True
        else:
            return False

    @classmethod
    def is_on_line(cls, pt1, pt2, pt3):
        """
        三点がa,b,c一直線に並ぶかどうか判定
        """
        vec_a = cls.vector(pt2, pt1)
        vec_b = cls.vector(pt3, pt1)
        return cls.is_parallel(vec_a, vec_b)

    @classmethod
    def split_triangle(cls, triangle, point):
        """
        三角形abcの内部の一点point:dから各頂点に直線を引いて三角形を分割する
        頂点上に乗る場合(分割なし),辺上(2分割),それ以外(3分割)
        """
        a = triangle[0]
        b = triangle[1]
        c = triangle[2]
        d = point
        # 頂点上
        if a[0:2] == d[0:2]:
            new_triangles = [[d, b, c]]
        elif b[0:2] == d[0:2]:
            new_triangles = [[a, d, c]]
        elif c[0:2] == d[0:2]:
            new_triangles = [[a, b, d]]
        # 辺上
        elif cls.is_on_line(a, b, d):
            new_triangles = [[d, b, c], [a, d, c]]
        elif cls.is_on_line(b, c, d):
            new_triangles = [[a, d, c], [a, b, d]]
        elif cls.is_on_line(c, a, d):
            new_triangles = [[d, b, c], [a, b, d]]
        else:
            new_triangles = [[d, b, c], [a, d, c], [a, b, d]]
        return new_triangles

    @classmethod
    def resplit_triangle(cls, tri1, tri2):
        """
        辺を共有する三角形(tri1, tri2)を分解して2種類の3角形の取り方を比較、
        安定な方を採用
        三角形を組み直した場合 ３要素目にFalseをreturnする
        """
        joint = [x for x in tri1 if x in tri2]
        not_joint = [x for x in tri1 if x not in tri2] + \
                        [x for x in tri2 if x not in tri1]
        judge = cls.which_lines(joint, not_joint)
        if joint == judge:
            return tri1, tri2, True
        alt1 = not_joint + [joint[0]]
        alt2 = not_joint + [joint[1]]

        return alt1, alt2, False

    @staticmethod
    def which_lines(line1, line2):
        """
        二線分の交点を算出
        安定な方をreturn
        """
        x1 = line1[0][0]
        y1 = line1[0][1]
        x2 = line1[1][0]
        y2 = line1[1][1]

        x3 = line2[0][0]
        y3 = line2[0][1]
        x4 = line2[1][0]
        y4 = line2[1][1]

        ksi = (y4-y3)*(x4-x1) - (x4-x3)*(y4-y1)
        eta = -(y2-y1)*(x4-x1)+(x2-x1)*(y4-y1)
        delta = (y4-y3)*(x2-x1)-(x4-x3)*(y2-y1)

        e1 = line1[0][2]
        e2 = line1[1][2]
        e3 = line2[0][2]
        e4 = line2[1][2]

        if delta**2 < 1e-16:  # 平行の場合
            return line1
        else:
            if 0 <= ksi/delta <= 1 and 0 <= eta/delta <= 1:  # 線分内にあるか
                energy_A = e1+(e2-e1)*ksi/delta
                energy_B = e4+(e3-e4)*eta/delta

                if energy_A == energy_B:
                    print('Both lines are same energy !!')
                    return line1

                elif energy_A > energy_B:
                    return line2

                elif energy_B > energy_A:
                    return line1
            return line1

        # 参考: 交点は以下で求められる
        # crosspoint = [x1+(x2-x1)*ksi/delta,
        #               y1+(y2-y1)*ksi/delta,
        #               energy_A, energy_B]


    @classmethod
    def is_inside(cls, triangle, point):
        """
        三角形triangle a,b,cの内部に点point dが位置するか判定
        辺上にある場合も内部と判定する
        """
        a, b, c = triangle
        d = point
        vec_b = cls.vector(b, a)
        vec_c = cls.vector(c, a)
        vec_d = cls.vector(d, a)
        coeff = cls.coeff_synthetic_vect(vec_b, vec_c, vec_d)
        judge = (0 <= coeff[0] and 0 <= coeff[1]
                 and 0 <= coeff[0]+coeff[1] <= 1)
        return judge

    @classmethod
    def separate_inout(cls, triangle, points):
        """三角形内部の点と外部の点に分別する"""
        inout = {True: [], False: []}
        for pt in points:
            inout[cls.is_inside(triangle, pt)].append(pt)
        return inout[True], inout[False]

    @classmethod
    def is_exist_more_stable(cls, triangle, points):
        """
        三角形内部の点のなかで三角energy平面よりもenergyが低い点がないか探す
        存在する場合を最安定の点をreturnする
        存在しない場合、Falseをreturn
        """
        in_points = cls.separate_inout(triangle, points)[0]
        energies = [[cls.relative_energy(triangle, x), x] for x in in_points]
        energies.sort()
        try:
            is_stable = energies[0][0] < 0
        except IndexError:
            return False
        out = {True: energies[0][1], False: False}
        return out[is_stable]

    @classmethod
    def local_stable(cls, triangle, points):
        """
        三角形のenergy基準から最もenergyが低い点をreturnする
        return: [relative_energy, point]
        """
        in_points = cls.separate_inout(triangle, points)[0]
        energies = [[cls.relative_energy(triangle, x), x] for x in in_points]
        energies.sort()
        return energies[0]

    @staticmethod
    def find_adjacent_triangles(src_tri, triangles):
        """
        src_triに対して辺を共有する三角形を探し出す
        新たに分割した三角形をsrc_triに取る場合、2つ以上はない (高々ひとつ)
        """
        out_tri = []
        for tri in triangles:
            judge = [x for x in tri if not x in src_tri]
            if len(judge) == 1:
                out_tri.append(tri)
        return out_tri

    @classmethod
    def most_stable(cls, base_tris, points):
        """
        local_stableを集めて最安定の点を選択
        """
        energies = [[cls.local_stable(x, points), x] for x in base_tris]
        energies.sort()
        relative_e = energies[0][0][0]
        point = energies[0][0][1]
        triangle = energies[0][1]
        return relative_e, point, triangle

    @classmethod
    def collect_base_triangles(cls, initial_triangle, original_points):
        """
        分割を繰り返し最安定の点を三点の一組で出力する
        点を追加する度に基底の三角形を取り直す必要がある
        """
        triangles = [initial_triangle]  # 分割した際に要素が増える
        if not original_points:
            return triangles
        bases = []
        points = copy.deepcopy(original_points)
        #  最安定を選択
        relative_e, new_stable, eq_tri = cls.most_stable(triangles, points)
        while relative_e < -1e-10:
            new_triangles = cls.split_triangle(eq_tri, new_stable)
            triangles.remove(eq_tri)
            # 新しい三角形に対してより安定な三角形の取り方がないか判定
            add_triangles = []
            for tri in new_triangles:
                adj_tri = cls.find_adjacent_triangles(tri, triangles)
                judge = True
                if adj_tri:
                    tri1, tri2, judge = cls.resplit_triangle(tri, adj_tri[0])
                if not judge:
                    triangles.remove(adj_tri[0])
                    triangles.append(tri2)
                    add_triangles.append(tri1)
                else:
                    add_triangles.append(tri)
            for tri in add_triangles:
                triangles.append(tri)
            relative_e, new_stable, eq_tri = cls.most_stable(triangles, points)
        return triangles

    @staticmethod
    def outline_bases(base_set):
        """
        basesをgnuplotで三角形が閉じるように出力
        """
        lines = ""
        for triangle in base_set:
            for points in triangle:
                lines += "\t".join([str(x) for x in points]) + "\n"
            lines += "\t".join(str(x) for x in triangle[0]) + "\n\n\n"
        return lines

    @classmethod
    def fromGround(cls, bases, points, withGS=False):
        """basesから構成されるenergy surfaceからのrelative energy"""
        pts = copy.deepcopy(points)
        from_gs = []
        for triangle in bases:
            in_points = cls.separate_inout(triangle, pts)[0]
            for point in in_points:
                relative = cls.relative_energy(triangle, point)
                from_gs.append({False: [point[0], point[1], relative],
                                True: [[point[0], point[1], relative],
                                       triangle]}[withGS])
                pts.remove(point)
        return from_gs

    @classmethod
    def get_ave_volume(cls, triangle, point):
        """
        triangleに対して組成点の平均volumeをreturn
        pointの4点目にvolumeを入れておく
        """
        # 3番目の要素に組成比を追加
        tri = copy.deepcopy(triangle)
        pt = copy.deepcopy(point)
        for tpt in tri:
            tpt.insert(2, 1 - sum(tpt[0:2]))
        pt.insert(2, 1 - sum(pt[0:2]))
        ave_v = sum([x[i]*pt[i]*x[4] for x in tri for i in range(3)])
        return ave_v

    @classmethod
    def ave_volume(cls, base, point):
        """
        与えられた3点(base)のvolumeに対するpoint組成での平均体積を算出
        format of base: [A1, A2, A3]=[[x1, y1, e1, v1], [x2, y2, e2, v2],...]
        """
        base_va = cls.vector(base[1], base[0])
        base_vb = cls.vector(base[2], base[0])
        vec_point = cls.vector(point, base[0])
        coeff = cls.coeff_synthetic_vect(base_va, base_vb, vec_point)
        ave_volume = coeff[0] * base_va[3] + coeff[1] * base_vb[3] + base[0][3]
        return ave_volume


if __name__ == '__main__':
    pass
