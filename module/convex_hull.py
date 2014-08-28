#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
convex_hullを作成
"""
import copy
import math
import pylab
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from commopy import Cabinet
from mpl_toolkits.mplot3d import proj3d
def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, a, b],
                     [0, 0, 0, zback]])


proj3d.persp_transformation = orthogonal_proj

def draw_convex_hull(ax, initial_base, not_bases, elements, zlim, **kwargs):
    """
    convex_hullを描く
    """
    bases = FindGS.collect_base_triangles(initial_base, not_bases)
    data = np.array(initial_base + not_bases)

    pt3d = PlotTriangularCoord(ax)

    pt3d.plt_triangle_dot(data, **kwargs)
    pt3d.outline_base_plot(bases)

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

    def make_flame(self, elements, zlim):
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
        self.ax.set_zlim(*zlim)

        self.ax.axis("off")
        #self.ax.text(0, 0, zlim[0]/2., 'Energy')
        self.ax.text(0.5, math.sqrt(3)/2+0.1, 0, elements[0], size=18)
        self.ax.text(-0.1, -0.05, 0, elements[1], size=18)
        self.ax.text(1.1, -0.05, 0, elements[2], size=18)

    def outline_base_plot(self, bases, **kwargs):
        """
        basesを三角形状に結んでプロット
        """
        for three_points in bases:
            triangle = np.array(three_points + [three_points[0]])
            x, y = alt_triangle(triangle[:, 0], triangle[:, 1])
            self.plt_line(x, y, triangle[:, 2], **kwargs)

    def plt_triangle_dot(self, data, **kwargs):
        x, y, z = (data[:, i] for i in range(0, 3))
        tri_x, tri_y = alt_triangle(x, y)
        self.ax.scatter3D(np.ravel(tri_x), np.ravel(tri_y), np.ravel(z),
                          **kwargs)



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
        format of base: [A1,A2,A3]=[[x1,y1,e1],[x2,y2,e2],...]
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
    def collect_base_triangles(cls, initial_triangle, original_points):
        """分割を繰り返し最安定の点を三点の一組で出力する"""
        triangles = [initial_triangle]  # 分割した際に要素が増える
        bases = []
        points = copy.deepcopy(original_points)
        while triangles:
            new_stable = cls.is_exist_more_stable(triangles[0], points)
            if not new_stable:
                in_points = cls.separate_inout(triangles[0], points)[0]
                bases.append(triangles[0])
                triangles.pop(0)
                for p in in_points:
                    points.remove(p)
            else:
                new_triangles = cls.split_triangle(triangles[0], new_stable)
                triangles.pop(0)
                for triangle in new_triangles:
                    triangles.append(triangle)
        return bases

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
    def fromGround(cls, bases, points):
        """basesから構成されるenergy surfaceからのrelative energy"""
        pts = copy.deepcopy(points)
        from_gs = []
        for triangle in bases:
            in_points = cls.separate_inout(triangle, pts)[0]
            for point in in_points:
                relative = cls.relative_energy(triangle, point)
                from_gs.append([point[0], point[1], relative])
                pts.remove(point)
        return from_gs


if __name__ == '__main__':
    pass
