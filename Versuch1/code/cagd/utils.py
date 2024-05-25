#!/usr/bin/python
from cagd.vec import Vec2

import numpy as np


# Solves the system of linear equations Ax = res
# where A is a tridiagonal matrix with diag2 representing the main diagonal
# diag1 and diag3 represent the lower and upper diagonal respectively
# All four parameters are vectors of size n
# The first element of diag1 and the last element of diag3 are ignored
# Therefore diag1[i], diag2[i] and diag3[i] are located on the same row of A
def solve_tridiagonal_equation(diag1, diag2, diag3, res):
    assert (len(diag1) == len(diag2) == len(diag3) == len(res))
    solution = None

    dim = len(res)

    diag1[0] = 0
    diag3[-1] = 0

    v = [0]
    y = [Vec2(0, 0)]

    for i in range(1, dim + 1):
        z = 1. / (diag2[i - 1] - diag1[i - 1] * v[i - 1])
        v.append(z * diag3[i - 1])
        y.append(z * (res[i - 1] - diag1[i - 1] * y[i - 1]))

    solution = [Vec2(0, 0)] * dim
    solution[-1] = y[-1]
    for i in range(dim - 2, -1, -1):
        solution[i] = y[i + 1] - v[i + 1] * solution[i + 1]

    return solution


# Solves the system of linear equations Ax = res
# where A is an almost tridiagonal matrix with diag2 representing the main diagonal
# diag1 and diag3 represent the lower and upper diagonal respectively
# All four parameters are vectors of size n
# The first element of diag1 and the last element of diag3 represent the top right and bottom left elements of A
# diag1[i], diag2[i] and diag3[i] are located on the same row of A
def solve_almost_tridiagonal_equation(diag1, diag2, diag3, res):
    assert (len(diag1) == len(diag2) == len(diag3) == len(res))

    n = len(diag1)

    # Initialize the vectors
    v = [0.0] * n
    y = [Vec2(0.0, 0.0)] * n
    s = [0.0] * n
    z = [0.0] * n

    # Forward Elimination
    z[0] = 1.0 / diag2[0]
    y[0] = res[0] * z[0]
    v[0] = -z[0] * diag3[0]
    s[0] = -diag1[0] * z[0]

    for k in range(1, n):
        z[k] = 1.0 / (diag2[k] + diag1[k] * v[k - 1])
        v[k] = -z[k] * diag3[k]
        y[k] = z[k] * (res[k] - diag1[k] * y[k - 1])
        s[k] = -diag1[k] * s[k - 1] * z[k]

    # Back substitution
    t = [0.0] * n
    w = [Vec2(0.0, 0.0)] * n

    t[-1] = 1.0
    w[-1] = Vec2(0.0, 0.0)

    for k in range(n-2, -1, -1):
        t[k] = v[k] * t[k+1] + s[k]
        w[k] = v[k] * w[k+1] + y[k]

    # Calculate xn
    a = diag1[-1] * t[-2] + diag3[-1] * t[0] + diag2[-1]
    xn_x = (res[-1].x - diag1[-1] * w[-2].x - diag3[-1] * w[0].x) / a
    xn_y = (res[-1].y - diag1[-1] * w[-2].y - diag3[-1] * w[0].y) / a
    xn = Vec2(xn_x, xn_y)

    solution = [Vec2(0.0, 0.0)] * n
    solution[-1] = Vec2(xn, xn)

    for k in range(n-1, -1, -1):
        solution[k] = Vec2(t[k] * xn.x + w[k].x, t[k] * xn.y + w[k].y)

    return solution