from cagd.utils import solve_almost_tridiagonal_equation
from cagd.vec import Vec2


def config1():
    diag1 = [9.4, 55.99, 1.58, 19.82, 19.09, 18.43, 17.37, 7.3]
    diag2 = [32.58, 87.73, 0.43, 74.03, 15.32, 63.97, 53.66, 42.14]
    diag3 = [10.37, 16.99, 36.54, 12.49, 35.86, 16.84, 31.36, 20.76]
    res = [Vec2(72.02, 4.25), Vec2(7.86, 2.45), Vec2(39.56, 43.8), Vec2(7.0, 6.44),
           Vec2(59.98, 25.17), Vec2(24.74, 73.71), Vec2(14.91, 67.16), Vec2(45.08, 10.13)]
    return diag1, diag2, diag3, res


def test_solve_almost_tridiagonal_equation():
    params = config1()
    points = solve_almost_tridiagonal_equation(*params)
    assert len(points) == 8
    assert 2.16 < points[0].x < 2.17
    assert 0.15 < points[0].y < 0.16
    assert 0.28 < points[1].x < 0.29
    assert -0.20 < points[1].y < -0.19
    assert -8.14 < points[2].x < -8.13
    assert 0.68 < points[2].y < 0.69
    assert 1.16 < points[3].x < 1.17
    assert 1.19 < points[3].y < 1.20
    assert 6.55 < points[4].x < 6.56
    assert -7.68 < points[4].y < -7.67
    assert -1.75 < points[5].x < -1.74
    assert 3.34 < points[5].y < 3.35
    assert 0.93 < points[6].x < 0.94
    assert 0.08 < points[6].y < 0.09
    assert -0.17 < points[7].x < -0.16
    assert 0.15 < points[7].y < 0.16
