from cagd.utils import solve_tridiagonal_equation
from cagd.vec import Vec2


def config1():
    diag1 = [9.4, 55.99, 1.58, 19.82, 19.09, 18.43, 17.37, 7.3]
    diag2 = [32.58, 87.73, 0.43, 74.03, 15.32, 63.97, 53.66, 42.14]
    diag3 = [10.37, 16.99, 36.54, 12.49, 35.86, 16.84, 31.36, 20.76]
    res = [Vec2(72.02, 4.25), Vec2(7.86, 2.45), Vec2(39.56, 43.8), Vec2(7.0, 6.44),
           Vec2(59.98, 25.17), Vec2(24.74, 73.71), Vec2(14.91, 67.16), Vec2(45.08, 10.13)]
    return diag1, diag2, diag3, res


def test_solve_tridiagonal_equation():
    params = config1()
    points = solve_tridiagonal_equation(*params)
    assert len(points) == 8
    assert 2.19 < points[0].x < 2.20
    assert 0.21 < points[0].y < 0.22
    assert 0.03 < points[1].x < 0.04
    assert -0.26 < points[1].y < -0.25
    assert -6.97 < points[2].x < -6.96
    assert 0.76 < points[2].y < 0.77
    assert 1.16 < points[3].x < 1.17
    assert 1.20 < points[3].y < 1.21
    assert 4.71 < points[4].x < 4.72
    assert -7.81 < points[4].y < -7.80
    assert -0.97 < points[5].x < -0.96
    assert 3.39 < points[5].y < 3.40
    assert -0.05 < points[6].x < -0.04
    assert 0.01 < points[6].y < 0.02
    assert 1.07 < points[7].x < 1.08
    assert 0.23 < points[7].y < 0.24
