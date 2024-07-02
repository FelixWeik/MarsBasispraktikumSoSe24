#!/usr/bin/python

from cagd.vec import Vec2, Vec3
from cagd.polyline import Polyline
import copy


class BezierCurve:
    def __init__(self, degree):
        assert degree >= 0
        self.degree = degree
        self.control_points = [None] * (degree + 1)
        self.color = "black"

    def set_control_point(self, index, val):
        assert 0 <= index <= self.degree
        self.control_points[index] = val

    def get_control_point(self, index):
        assert 0 <= index <= self.degree
        return self.control_points[index]

    # Evaluates the curve at t
    def evaluate(self, t):
        return self._de_casteljau(t, 1)[0]

    # Evaluates tangent at t
    def tangent(self, t):
        last_two_ctrl_pts = self._de_casteljau(t, 2)
        a = last_two_ctrl_pts[0]
        b = last_two_ctrl_pts[1]
        return b - a

    # calculates the normal at t
    def normal(self, t):
        pass

    # Syntactic sugar so bezier curve can be evaluated as curve(t)
    # instead of curve.evaluate(t)
    def __call__(self, t):
        return self.evaluate(t)

    # Calculates the de Casteljau scheme until the column only has stop elements
    def _de_casteljau(self, t, stop):
        assert (stop >= 1)
        column = self.control_points
        while len(column) > stop:
            new_column = [None for i in range(len(column) - 1)]
            for i in range(len(new_column)):
                new_column[i] = (1 - t) * column[i] + t * column[i + 1]
            column = new_column
        return column

    def get_color(self):
        return self.color

    def set_color(self, color):
        self.color = color

    # Calculates the bezier representation of the derivative
    def get_derivative(self):
        pass

    def get_axis_aligned_bounding_box(self):
        min_vec = copy.copy(self.control_points[0])
        max_vec = copy.copy(self.control_points[0])
        for p in self.control_points:
            if p.x < min_vec.x:
                min_vec.x = p.x
            if p.y < min_vec.y:
                min_vec.y = p.y
            if p.x > max_vec.x:
                max_vec.x = p.x
            if p.y > max_vec.y:
                max_vec.y = p.y
        return min_vec, max_vec

    def draw(self, scene, num_samples):
        p0 = self(0)
        for i in range(1, num_samples + 1):
            t = i / num_samples
            p1 = self(t)
            scene.draw_line(p0, p1, self.color)
            p0 = p1

    def get_polyline_from_control_points(self):
        pl = Polyline()
        for p in self.control_points:
            pl.append_point(p)
        return pl


class BezierSurface:
    # Creates a bezier surface of degrees n,m
    # The degree parameter is a tuple (n,m)
    def __init__(self, degree):
        d1, d2 = degree
        assert (d1 >= 0 and d2 >= 0)
        self.degree = degree
        self.control_points = [[None for i in range(d2 + 1)] for j in range(d1 + 1)]
        white = (1, 1, 1)
        self.color = (white, white, white, white)
        self.curvature = (None, None, None, None)

    def set_control_point(self, index1, index2, val):
        assert 0 <= index1 <= self.degree[0]
        assert 0 <= index2 <= self.degree[1]
        self.control_points[index1][index2] = val

    def get_control_point(self, index1, index2):
        assert 0 <= index1 <= self.degree[0]
        assert 0 <= index2 <= self.degree[1]
        return self.control_points[index1][index2]

    def evaluate(self, t1, t2):
        return self._de_casteljau(t1, t2, (1, 1))[0][0]

    # Sets the colors at the corners
    # c00 is the color at u=v=0, c01 is the color at u=0 v=1, etc.
    # A color is a tuple (r,g,b) with values between 0 an 1
    def set_colors(self, c00, c01, c10, c11):
        self.color = (c00, c01, c10, c11)

    # Sets the curvature at the corners
    # c00 is the curvature at u=v=0, c01 is the curvature at u=0 v=1, etc
    def set_curvature(self, c00, c01, c10, c11):
        self.curvature = (c00, c01, c10, c11)

    def __call__(self, t):
        t1, t2 = t
        return self.evaluate(t1, t2)

    def _de_casteljau(self, t1, t2, stop):
        s1, s2 = stop
        d1, d2 = self.degree
        assert (s1 >= 1 and s2 >= 1)
        d1 += 1  # number of control points in each direction
        d2 += 1

        # Apply the de Casteljau scheme in one direction,
        # ie, reduce dimension from (d1, d2) to (s1, d2)
        column = self.control_points
        while d1 > s1:
            d1 -= 1
            new_column = [[None for i in range(d2)] for j in range(d1)]
            for i in range(d1):
                for j in range(d2):
                    new_column[i][j] = (1 - t1) * column[i][j] + t1 * column[i + 1][j]
            column = new_column

        # Apply the de Casteljau scheme in the other direction,
        # ie, reduce dimension from (s1, d2) to (s1, s2)
        while d2 > s2:
            d2 -= 1
            new_column = [[None for i in range(d2)] for j in range(d1)]
            for i in range(d1):
                for j in range(d2):
                    new_column[i][j] = (1 - t2) * column[i][j] + t2 * column[i][j + 1]
            column = new_column

        return column

    def normal(self, t1, t2):
        pass

    def get_derivative(self, direction):
        pass

    def subdivide(self, t1, t2):
        b0, b1 = self._subdivide_u(t1)
        b00, b01 = b0._subdivide_v(t2)
        b10, b11 = b1._subdivide_v(t2)
        return [b00, b01, b10, b11]

    def _subdivide_u(self, t):
        du, dv = self.degree
        left = BezierSurface((du, dv))
        right = BezierSurface((du, dv))
        for k in range(du + 1):
            pts = self._de_casteljau(t, 0, (du - k + 1, dv + 1))
            left.control_points[k] = pts[0]
            right.control_points[-(k + 1)] = pts[-1]
        return left, right

    def _subdivide_v(self, t):
        du, dv = self.degree
        left = BezierSurface((du, dv))
        right = BezierSurface((du, dv))
        for k in range(dv + 1):
            pts = self._de_casteljau(0, t, (du + 1, dv - k + 1))
            for i in range(du + 1):
                left.control_points[i][k] = pts[i][0]
                right.control_points[i][-(k + 1)] = pts[i][-1]
        return left, right


class BezierPatches:
    CURVATURE_GAUSSIAN = 0
    CURVATURE_AVERAGE = 1
    CURVATURE_PRINCIPAL_MAX = 2  # Maximale Hauptkruemmung
    CURVATURE_PRINCIPAL_MIN = 3  # Minimale Hauptkruemmung
    COLOR_MAP_LINEAR = 4
    COLOR_MAP_CUT = 5
    COLOR_MAP_CLASSIFICATION = 6

    def __init__(self):
        self.patches = []

    def __len__(self):
        return len(self.patches)

    def __getitem__(self, p):
        return self.patches[p]

    def __setitem__(self, i, p):
        self.patches[i] = p

    def __delitem__(self, p):
        del self.patches[p]

    def __iter__(self):
        return iter(self.patches)

    def append(self, p):
        self.patches.append(p)

    # Refines patches by subdividing each patch into four new patches
    # There are 4^num times more patches after calling this function
    def refine(self, num):
        for i in range(num):
            new_patches = BezierPatches()
            for p in self:
                new = p.subdivide(0.5, 0.5)
                for n in new:
                    new_patches.append(n)
            self.patches = new_patches

   def visualize_curvature(self, curvature_mode, color_map):
        # coordinates of each corner
        corners = [(0, 0), (0, 1), (1, 0), (1, 1)]

        # set curvature to maximum or minimum values to update afterwards
        max_curve_value = float('-inf')
        min_curve_value = float('inf')

        # for each patch
        for patch in self.patches:
            # get the degree u and v of the patch
            u = patch.degree[0]
            v = patch.degree[1]

            # 1st partial derivatives
            du = BezierSurface([u - 1, v])
            dv = BezierSurface([u, v - 1])

            # 2nd partial derivatives
            duu = BezierSurface([u - 2, v])
            dvv = BezierSurface([u, v - 2])
            dvu = BezierSurface([u - 1, v - 1])

            # set control points for the 1st partial derivatives
            for i in range(u):
                for j in range(v + 1):
                    # rate of change between consecutive control points in the u direction, scaled by the degree u
                    du.control_points[i][j] = (patch.control_points[i + 1][j] - patch.control_points[i][j]) * u
                    # rate of change between consecutive control points in the v direction, scaled by the degree v
                    dv.control_points[j][i] = (patch.control_points[j][i + 1] - patch.control_points[j][i]) * v

            # set control points for the 2nd partial derivatives
            # partial derivative of du with respect to u
            for i in range(u - 1):
                for j in range(v + 1):
                    # rate of change between consecutive control points in the u direction, scaled by the degree u - 1
                    duu.control_points[i][j] = (du.control_points[i + 1][j] - du.control_points[i][j]) * (u - 1)
            # partial derivative of dv with respect to v
            for i in range(u):
                for j in range(v - 1):
                    # rate of change between consecutive control points in the v direction, scaled by the degree v - 1
                    dvv.control_points[j][i] = (dv.control_points[j][i + 1] - dv.control_points[j][i]) * (v - 1)
            # mixed partial derivative of dv with respect to u
            for i in range(u - 1):
                for j in range(v):
                    # rate of change between consecutive control points in the u direction, scaled by the degree u - 1
                    dvu.control_points[i][j] = (dv.control_points[i + 1][j] - dv.control_points[i][j]) * (u - 1)

            curvatures = [0, 0, 0, 0]
            # for each corner
            for corner_index, (c1, c2) in enumerate(corners):
                # get the control points of each partial derivative of the corner
                corner_du = du.control_points[c1 * du.degree[0]][c2 * du.degree[1]]
                corner_dv = dv.control_points[c1 * dv.degree[0]][c2 * dv.degree[1]]
                corner_duu = duu.control_points[c1 * duu.degree[0]][c2 * duu.degree[1]]
                corner_dvv = dvv.control_points[c1 * dvv.degree[0]][c2 * dvv.degree[1]]
                corner_dvu = dvu.control_points[c1 * dvu.degree[0]][c2 * dvu.degree[1]]

                # bu.bu
                E = corner_du.dot(corner_du)
                # bu.bv
                F = corner_du.dot(corner_dv)
                # bv.bv
                G = corner_dv.dot(corner_dv)

                # calculate the normal vector
                N = corner_du.cross(corner_dv).div(corner_du.cross(corner_dv).length())

                # normal.buu
                e = N.dot(corner_duu)
                # normal.buv
                f = N.dot(corner_dvu)
                # normal.bvv
                g = N.dot(corner_dvv)

                # calculate the different curvature methods
                gaussian = (e * g - f * f) / (E * G - F * F)

                average = 0.5 * (e * G - 2 * f * F + g * E) / (E * G - F * F)

                principal1 = average + math.sqrt(average * average - gaussian)
                principal2 = average - math.sqrt(average * average - gaussian)

                # choose the curvature method according to the parameter
                if curvature_mode == self.CURVATURE_GAUSSIAN:
                    curvatures[corner_index] = gaussian
                elif curvature_mode == self.CURVATURE_AVERAGE:
                    curvatures[corner_index] = average
                elif curvature_mode == self.CURVATURE_PRINCIPAL_MAX:
                    curvatures[corner_index] = max(principal1, principal2)
                elif curvature_mode == self.CURVATURE_PRINCIPAL_MIN:
                    curvatures[corner_index] = min(principal1, principal2)

                # update the maximum and minimum curvature values
                max_curve_value = max(max_curve_value, curvatures[corner_index])
                min_curve_value = min(min_curve_value, curvatures[corner_index])

            # set the curvatures to the patch
            patch.set_curvature(*curvatures)

        # for each patch
        for patch in self.patches:
            # a color is of the form (r, g, b) with values between 0 an 1
            colors = [(0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)]

            # for each corner
            for corner_index, (c1, c2) in enumerate(corners):
                # get the curvature value
                corner_curve_value = patch.curvature[corner_index]

                # choose the color map method according to the parameter
                if color_map == self.COLOR_MAP_CUT:
                    # if the curvature is less than 0, the color is set to 0.
                    # if it is more than 1, the color is set to 1.
                    color_curve_value = min(1, max(0, corner_curve_value))
                elif color_map == self.COLOR_MAP_LINEAR:
                    # the minimum value of curvature will be mapped to a 0
                    # the maximum value of curvature will be mapped to a 1
                    # all the values in between are set linearly
                    color_curve_value = (corner_curve_value - min_curve_value) / (max_curve_value - min_curve_value)
                elif color_map == self.COLOR_MAP_CLASSIFICATION:
                    # this classifies the curvature value into three categories:
                    # 0 for negative curvature values (less than -0.0001).
                    # 0.5 for near-zero curvature values (between -0.0001 and 0.0001).
                    # 1 for positive curvature values (greater than 0.0001)
                    color_curve_value = 0 if corner_curve_value < -0.0001 else (1 if corner_curve_value > 0.0001 else 0.5)

                colors[corner_index] = BezierPatches.get_rgb_from_color_value(color_curve_value)

            # update the color of the patch
            patch.set_colors(*colors)

    def get_rgb_from_color_value(value):
        # from the value, we return a tuple lile (r, g, b)
        if 0 <= value <= 0.25:
            # for values between 0 and 0.25, the color transitions from blue to cyan
            return 0, 4 * value, 1
        elif 0.25 < value <= 0.5:
            # for values between 0.25 and 0.5, the color transitions from cyan to green
            return 0, 1, 2 - 4 * value
        elif 0.5 < value <= 0.75:
            # for values between 0.5 and 0.75, the color transitions from green to yellow
            return 4 * value - 2, 1, 0
        elif 0.75 < value <= 1:
            # for values between 0.75 and 1, the color transitions from yellow to red
            return 1, 4 - 4 * value, 0

    def f1(self, curvature_value):
        normalized_curvature = (curvature_value - min(curvature_value.values())) / (max(curvature_value.values()) - min(curvature_value.values()))
        r = int(255 * normalized_curvature)
        b = int(255 * (1 - normalized_curvature))
        color = (r, 0, b)
        return color

def get_derivative(self, direction):
        if direction not in ['u', 'v']:
            raise ValueError("Invalid direction. Use 'u' or 'v'.")

        n = len(self.control_points) - 1 
        derivative_control_points = []

        if direction == 'u':
            for i in range(n):
                row = []
                for j in range(len(self.control_points[i])):
                    du = n * (self.control_points[i + 1][j] - self.control_points[i][j])
                    row.append(du)
                derivative_control_points.append(row)
        elif direction == 'v':
            for i in range(len(self.control_points)):
                row = []
                for j in range(n):
                    dv = n * (self.control_points[i][j + 1] - self.control_points[i][j])
                    row.append(dv)
                derivative_control_points.append(row)

        return derivative_control_points

        def f2(self, curvature_value):
        pass

    def f3(self, curvature_value):
        pass

    def export_off(self):
        def export_point(p):
            return str(p.x) + " " + str(p.y) + " " + str(p.z)

        def export_colors(c):
            s = ""
            for x in c:
                s += str(x)
                s += " "
            s += "1"  # opacity
            return s

        s = "CBEZ333\n"
        for patch in self:
            # coordinates
            for row in patch.control_points:
                for point in row:
                    s += export_point(point)
                    s += "\n"

            # colors
            s += export_colors(patch.color[0])
            s += "\n"
            s += export_colors(patch.color[2])
            s += "\n"
            s += export_colors(patch.color[1])
            s += "\n"
            s += export_colors(patch.color[3])
            s += "\n"
            s += "\n"

        return s

    def export_standard_off(self):
        def export_numbers(degree):
            d1, d2 = degree
            n_v = (d1 + 1) * (d2 + 1)
            n_f = d1 * d2
            n_e = d1 * (d2 + 1) + d2 * (d1 + 1)
            patch_num = len(self)
            return str(n_v * patch_num) + " " + str(n_f * patch_num) + " " + str(n_e * patch_num) + "\n"

        def export_vertex(v):
            return str(v.x) + " " + str(v.y) + " " + str(v.z) + "\n"

        def export_face(f_vs):
            s = "4  "
            for f_v in f_vs:
                s += str(f_v)
                s += " "
            return s

        def avg_color(f_cs):
            avg_f_c = [0, 0, 0]
            for f_c in f_cs:
                avg_f_c[0] += f_c[0]
                avg_f_c[1] += f_c[1]
                avg_f_c[2] += f_c[2]
            return (round(255 * avg_f_c[0] / 4),
                    round(255 * avg_f_c[1] / 4),
                    round(255 * avg_f_c[2] / 4))

        def export_color(f_c):
            return " " + str(f_c[0]) + " " + str(f_c[1]) + " " + str(f_c[2]) + "\n"

        s = "OFF\n\n"
        s += export_numbers(self[0].degree)

        for patch in self:
            for v_row in patch.control_points:
                for vertex in v_row:
                    s += export_vertex(vertex)
        start_v = 0
        for patch in self:
            cps = patch.control_points
            row_num = len(cps)
            col_num = len(cps[0])
            for row in range(row_num - 1):
                for col in range(col_num - 1):
                    first_v = start_v + row * col_num + col
                    f_vertices = (first_v, first_v + 1, first_v + 1 + col_num, first_v + col_num)
                    s += export_face(f_vertices)
                    s += export_color(avg_color(patch.color))
            start_v += 16
        return s
