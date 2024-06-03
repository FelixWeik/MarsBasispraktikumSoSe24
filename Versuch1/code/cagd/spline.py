#! /usr/bin/python

import math

import numpy as np
from cagd.vec import Vec2, Vec3
from cagd.polyline import Polyline
from cagd.bezier import BezierSurface, BezierPatches
import cagd.utils as utils
import copy

class Spline:
    # Interpolation modes
    INTERPOLATION_GIVEN_KNOTS = 0
    INTERPOLATION_EQUIDISTANT = 1
    INTERPOLATION_CHORDAL = 2
    INTERPOLATION_CENTRIPETAL = 3
    INTERPOLATION_FOLEY = 4

    def __init__(self, degree):
        assert (degree >= 1)
        self.degree = degree
        self.periodic = False
        self.knots = None
        self.control_points = []
        self.color = "black"

    # Checks if the number of knots, control points and degree define a valid spline
    def validate(self):
        knots = self.knots.validate()
        points = len(self.knots) == len(self.control_points) + self.degree + 1
        return knots and points

    def evaluate(self, t):
        a, b = self.support()
        assert (a <= t <= b)
        if t == self.knots[len(self.knots) - self.degree - 1]:
            # the spline is only defined on the interval [a, b)
            # it is useful to define self(b) as lim t->b self(t)
            t = t - 0.000001
        return self.de_boor(t, 1)[0]

    # Returns the interval [a, b) on which the spline is supported
    def support(self):
        return self.knots[self.degree], self.knots[len(self.knots) - self.degree - 1]

    def __call__(self, t):
        return self.evaluate(t)

    def tangent(self, t):
        # same as evaluate() method #############################
        a, b = self.support()
        assert (a <= t <= b)
        if t == self.knots[len(self.knots) - self.degree - 1]:
            # the spline is only defined on the interval [a, b)
            # it is useful to define self(b) as lim t->b self(t)
            t = t - 0.000001
        # ends here #############################################

        # apply De Boor algorithm, to evaluate B-splines.
        # return a list of 2 points that influence the spline at t.
        pts = self.de_boor(t, 2)

        # computes the difference between two points,
        # giving a vector 'tangent_direction' that represents the direction of the tangent.
        tangent_direction = pts[1] - pts[0]

        # compute the length by using the dot product
        length = math.sqrt(tangent_direction.dot(tangent_direction))

        # we normalize the direction vector by dividing its components by its length.
        return Vec2(tangent_direction.x / length, tangent_direction.y / length)

    def get_color(self):
        return self.color

    def set_color(self, color):
        self.color = color

    # Calculates the de_boor scheme at a given value t
    # Stops when the column is only "stop" elements long
    # Returns that column as a list
    def de_boor(self, t, stop):
        index = self.knots.knot_index(t)
        points = []

        if (index - self.degree) < 0 or index + 1 > len(self.control_points):
            raise AssertionError("Index out of range")

        for i in range(index - self.degree, index + 1):
            points.append(self.control_points[i])

        knot_vector = []
        if (index - self.degree + 1) < 0 or index + self.degree + 1 > len(self.knots):
            raise AssertionError("Index out of range")

        for i in range(index - self.degree + 1, index + self.degree + 1):
            knot_vector.append(self.knots[i])
        k = 1
        while points.__len__() > 1 and points.__len__() > stop:

            next_points = []
            for j in range(0, points.__len__() - 1):
                alpha = (t - knot_vector[j + k - 1]) / (
                            knot_vector[j + self.degree] - knot_vector[j + k - 1])
                next_points.append((1 - alpha) * points[j] + alpha * points[j + 1])

            points = next_points
            k += 1
        return points

    # Adjusts the control points such that it represents the same function,
    # but with an added knot
   def insert_knot(self, t):
        # find the appropriate index in the knot vector where the new knot should be inserted.
        index = self.knots.knot_index(t)

        # compute new control points
        new_control_points = self.de_boor(t, self.degree)

        # determine the 1st point to be replaced
        first_point = index - self.degree + 1

        self.control_points = (
            self.control_points[:first_point] # keep the initial control points
            + new_control_points # add the new control points where they should be
            + self.control_points[index:] # then add the remaining intial control points
        )

        # insert the node in sorted order.
        self.knots.insert(t)

    def get_axis_aligned_bounding_box(self):
        min_vec = copy.copy(self.control_points[0])
        max_vec = copy.copy(self.control_points[0])
        for p in self.control_points:
            # print("comparing {0} to {1} and {2}".format(p, min_vec, max_vec))
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
        i = self.degree - 1
        while i < len(self.knots) - self.degree - 2:
            i += 1
            k0 = self.knots[i]
            k1 = self.knots[i + 1]
            if k0 == k1:
                continue
            p0 = self(k0)
            for j in range(1, num_samples + 1):
                t = k0 + j / num_samples * (k1 - k0)
                p1 = self(t)
                scene.draw_line(p0, p1, self.color)
                p0 = p1

    def get_polyline_from_control_points(self):
        pl = Polyline()
        for p in self.control_points:
            pl.append_point(p)
        return pl

    # Generates a spline that interpolates the given points using the given mode
    # kts is only used as given knots in the mode: INTERPOLATION_GIVEN_KNOTS
    # Returns that spline object
    @classmethod
    def interpolate_cubic(cls, mode, points, kts=None):
        spline = Spline(3)
        match mode:
            case cls.INTERPOLATION_GIVEN_KNOTS:
                spline.knots = copy.deepcopy(kts)

            case cls.INTERPOLATION_EQUIDISTANT:
                knots = Knots(len(points) + 6)
                for i in range(knots.knots.__len__()):
                    knots[i] = np.clip(i - 3, 0, len(points) - 1)
                spline.knots = knots
            case cls.INTERPOLATION_CHORDAL:
                knots = Knots(len(points) + 6)
                for i in range(4):
                    knots[i] = 0

                for i in range(len(points) - 1):
                    diff = points[i + 1] - points[i]
                    knots[i + 4] = knots[i + 3] + math.sqrt(diff.dot(diff))

                for i in range(len(points) + 3, len(points) + 6):
                    knots[i] = knots[i - 1]
                spline.knots = knots
            case cls.INTERPOLATION_CENTRIPETAL:
                knots = Knots(len(points) + 6)
                for i in range(4):
                    knots[i] = 0

                for i in range(len(points) - 1):
                    diff = points[i + 1] - points[i]
                    knots[i + 4] = knots[i + 3] + math.sqrt(math.sqrt(diff.dot(diff)))

                for i in range(len(points) + 3, len(points) + 6):
                    knots[i] = knots[i - 1]
                spline.knots = knots

        dim = len(points) + 2

        alpha = []
        beta = []
        gamma = []
        for i in range(2, dim - 2):
            alpha.append((spline.knots[i + 2] - spline.knots[i]) / (spline.knots[i + 3] - spline.knots[i]))
            beta.append((spline.knots[i + 2] - spline.knots[i + 1]) / (spline.knots[i + 3] - spline.knots[i + 1]))
            gamma.append((spline.knots[i + 2] - spline.knots[i + 1]) / (spline.knots[i + 4] - spline.knots[i + 1]))

        col1 = [0] * dim
        col1[0] = 0
        col1[1] = -1
        col1[-1] = 0
        col1[-2] = -1 + gamma[-1]

        col2 = [0] * dim
        col2[0] = 1
        col2[1] = 1 + alpha[0]
        col2[-1] = 1
        col2[-2] = -gamma[-1] + 2

        col3 = [0] * dim
        col3[0] = 0
        col3[1] = -alpha[0]
        col3[-1] = 0
        col3[-2] = -1

        for i in range(2, dim - 2):
            col1[i] = (1 - beta[i - 2]) * (1 - alpha[i - 2])
            col2[i] = (1 - beta[i - 2]) * alpha[i - 2] + beta[i - 2] * (1 - gamma[i - 2])
            col3[i] = beta[i - 2] * gamma[i - 2]

        res = [Vec2(0, 0)] * dim
        res[0] = points[0]
        res[-1] = points[-1]

        for i in range(2, dim - 2):
            res[i] = points[i - 1]

        spline.control_points = utils.solve_tridiagonal_equation(col1, col2, col3, res)

        return spline



    # Generates a spline that interpolates the given points and fulfills the definition
    # of a periodic spline with equidistant knots
    # Returns that spline object
    @classmethod
    def interpolate_cubic_periodic(cls, points):
        degree = 3
        spline = Spline(degree)
        spline.periodic = True

        # Create a periodic knot vector with equidistant knots
        knots = Knots(12)
        for i in range(len(knots)):
            knots[i] = i
        spline.knots = knots

        # Setting up the coefficients for the almost tridiagonal system
        alpha = 1.0 / 6.0
        beta = 4.0 / 6.0
        gamma = 1.0 / 6.0

        dim = len(points)

        diag1 = [gamma] * dim
        diag2 = [beta] * dim
        diag3 = [alpha] * dim

        res = [Vec2(0, 0)] * dim
        for i in range(dim):
            res[i] = points[i]

        # Solving the almost tridiagonal system to get control points
        spline.control_points = utils.solve_almost_tridiagonal_equation(diag1, diag2, diag3, res)
        for i in range(degree):
            spline.control_points.append(spline.control_points[i])

        return spline

    # For splines of degree 3, generate a parallel spline with distance dist
    # The returned spline is off from the exact parallel by at most eps
    def generate_parallel(self, dist, eps):
        assert (self.degree == 3)
        if dist == 0:
            return self

        para_spline = None
        return para_spline

    # Generates a rotational surface by rotating the spline around the z axis
    # the spline is assumed to be on the xz-plane
    # num_samples refers to the number of interpolation points in the rotational direction
    # Returns a spline surface object in three dimensions
    def generate_rotation_surface(self, num_samples):
        pass


class SplineSurface:
    # The two directions of the parameter space
    DIR_U = 0
    DIR_V = 1

    # Creates a spline of degrees n,m
    # Degree is a tuple (n,m)
    def __init__(self, degree):
        du, dv = degree
        assert (du >= 1 and dv >= 1)
        self.degree = degree
        self.periodic = (False, False)
        self.knots = (None, None)  # tuple of both knot vectors
        self.control_points = [[]]  # 2dim array of control points

    # Checks if the number of knots, control points and degree define a valid spline
    def validate(self):
        if len(self.control_points) == 0:
            return False
        k1, k2 = self.knots
        d1, d2 = self.degree
        knots12 = k1.validate() and k2.validate()
        p1 = len(self.control_points)
        p2 = len(self.control_points[0])
        points1 = len(k1) == p1 + d1 + 1
        points2 = len(k2) == p2 + d2 + 1
        return knots12 and points1 and points2

    def evaluate(self, u, v):
        s1, s2 = self.support()
        a, b = s1
        c, d = s2
        assert (a <= u <= b and c <= v <= d)
        if u == b:
            u = u - 0.000001
        if v == d:
            v = v - 0.000001
        t = (u, v)
        return self.de_boor(t, (1, 1))[0][0]

    # Return nested tuple ((a,b), (c,d))
    # The spline is supported in (u,v) \in [a,b)x[c,d]
    def support(self):
        k1, k2 = self.knots
        d1, d2 = self.degree
        s1 = (k1[d1], k1[len(k1) - d1 - 1])
        s2 = (k2[d2], k2[len(k2) - d2 - 1])
        return s1, s2

    def __call__(self, u, v):
        return self.evaluate(u, v)

    # Calculates the de boor scheme at t = (u,v)
    # until there are only stop = (s1, s2) elements left
    def de_boor(self, t, stop):
        d1, d2 = self.degree
        k1, k2 = self.knots
        s1, s2 = stop
        u, v = t
        m1 = len(self.control_points)
        m2 = len(self.control_points[0])

        new_rows = [None] * m1
        for row in range(m1):
            spl = Spline(d2)
            spl.knots = k2
            spl.control_points = self.control_points[row]
            new_rows[row] = spl.de_boor(v, s2)

        new_pts = [None] * s2
        for col in range(s2):
            spl = Spline(d1)
            spl.knots = k1
            ctrl_pts = [new_rows[i][col] for i in range(m1)]
            spl.control_points = ctrl_pts
            new_pts[col] = spl.de_boor(u, s1)

        return new_pts

    def insert_knot(self, direction, t):
        if direction == self.DIR_U:
            self._insert_knot_u(t)
        elif direction == self.DIR_V:
            self._insert_knot_v(t)
        else:
            assert False

    def _insert_knot_v(self, t):
        du, dv = self.degree
        pu, pv = self.periodic
        ku, kv = self.knots
        nu = len(self.control_points)
        nv = len(self.control_points[0])
        for i in range(nu):
            row = self.control_points[i]
            spl = Spline(dv)
            spl.control_points = copy.copy(row)
            spl.knots = copy.deepcopy(kv)
            spl.periodic = pv
            spl.insert_knot(t)
            self.control_points[i] = spl.control_points
            self.knots = (ku, spl.knots)

    def _insert_knot_u(self, t):
        du, dv = self.degree
        pu, pv = self.periodic
        ku, kv = self.knots
        nu = len(self.control_points)
        nv = len(self.control_points[0])
        new_control_points = [[None for i in range(nv)] for j in range(nu + 1)]
        for i in range(nv):
            col = [self.control_points[j][i] for j in range(nu)]
            spl = Spline(du)
            spl.control_points = col
            spl.knots = copy.deepcopy(ku)
            spl.periodic = pu
            spl.insert_knot(t)
            for j in range(nu + 1):
                new_control_points[j][i] = spl.control_points[j]
            self.knots = (spl.knots, kv)
        self.control_points = new_control_points

    # Build bezier patches based on the spline with multiple knots
    # and control points sitting also as bezier points.
    def to_bezier_patches(self):
        patches = BezierPatches()
        return patches


class Knots:
    # Creates a knots array with n elements
    def __init__(self, n):
        self.knots = [None] * n

    def validate(self):
        prev = None
        for k in self.knots:
            if k is None:
                return False
            if prev is not None:
                if k < prev:
                    return False
            prev = k
        return True

    def __len__(self):
        return len(self.knots)

    def __getitem__(self, i):
        return self.knots[i]

    def __setitem__(self, i, v):
        self.knots[i] = v

    def __delitem__(self, i):
        del self.knots[i]

    def __iter__(self):
        return iter(self.knots)

    def insert(self, t):
        i = 0
        while self[i] < t:
            i += 1
        self.knots.insert(i, t)

    def knot_index(self, v):

        if self.knots is None or self.knots.__len__() == 0:
            return None;

        last = self.knots[-1]
        if last < v:
            return None

        i = 0

        if last == v:
            while i < self.__len__() and self.knots[i] < v:
                i += 1
            if i > 0:
                return i - 1
            else:
                return None

        while i < self.__len__() and self.knots[i] <= v:
            i += 1
        if i < 1:
            return None
        return i - 1
