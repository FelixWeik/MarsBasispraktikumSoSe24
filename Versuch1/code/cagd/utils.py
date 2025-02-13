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
    solution = None

    dim = len(res)

    v = [0] * dim
    y = [Vec2(0, 0)] * dim
    s = [0] * dim

    s[0] = 1

    for i in range(1, dim):
        z = 1 / (diag2[i - 1] + diag1[i - 1] * v[i - 1])
        v[i] = -z * diag3[i - 1]
        y[i] = z * (res[i - 1] - diag1[i - 1] * y[i - 1])
        s[i] = -diag1[i - 1] * s[i - 1] * z

    t = [0] * dim
    w = [Vec2(0, 0)] * dim
    t[-1] = 1

    for i in range(dim - 2, -1, -1):
        t[i] = v[i + 1] * t[i + 1] + s[i + 1]
        w[i] = v[i + 1] * w[i + 1] + y[i + 1]

    x_n = (res[-1] - diag3[-1] * w[0] - diag1[-1] * w[-2]) * (1.0 / (diag3[-1] * t[0] + diag1[-1] * t[-2] + diag2[-1]))

    x = [Vec2(0, 0)] * dim
    x[-1] = x_n

    for i in range(dim - 2, -1, -1):
        x[i] = t[i] * x_n + w[i]

    return x
