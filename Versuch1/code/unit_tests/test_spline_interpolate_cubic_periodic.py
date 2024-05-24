from cagd.spline import Spline
from cagd.vec import Vec2


def config1():
    pts = [Vec2(0.2, 3), Vec2(5, 1), Vec2(3, -0.2), Vec2(1, 2), Vec2(1, 3)]
    return pts


def test_interpolate_periodic_degree():
    points = config1()
    spl = Spline.interpolate_cubic_periodic(points)
    assert spl.degree == 3, "Incorrect degree."


def test_interpolate_periodic_knot_amount():
    points = config1()
    spl = Spline.interpolate_cubic_periodic(points)
    kts = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    assert len(spl.knots) == len(kts), "Incorrect number of Knots."


def test_interpolate_periodic_knot_values():
    points = config1()
    spl = Spline.interpolate_cubic_periodic(points)
    kts = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    assert len(spl.knots) == len(kts), "Incorrect number of Knots."
    for n, (i, j) in enumerate(zip(spl.knots, kts)):
        assert i == j, f"Wrong knot value at pos {n}."


def test_interpolate_periodic_control_points_amount():
    points = config1()
    spl = Spline.interpolate_cubic_periodic(points)
    assert len(spl.control_points) == 8, "Incorrect number of control points."


def test_interpolate_periodic_control_points_values():
    points = config1()
    spl = Spline.interpolate_cubic_periodic(points)
    assert -2.02 < spl.control_points[0].x < -2.01, "Incorrect x coordinate of control point at position 0."
    assert 3.52 < spl.control_points[0].y < 3.53, "Incorrect y coordinate of control point at position 0."
    assert 7.36 < spl.control_points[1].x < 7.37, "Incorrect x coordinate of control point at position 1."
    assert 0.90 < spl.control_points[1].y < 0.91, "Incorrect y coordinate of control point at position 1."
    assert 2.56 < spl.control_points[2].x < 2.57, "Incorrect x coordinate of control point at position 2."
    assert -1.17 < spl.control_points[2].y < -1.16, "Incorrect y coordinate of control point at position 2."
    assert 0.38 < spl.control_points[3].x < 0.39, "Incorrect x coordinate of control point at position 3."
    assert 2.54 < spl.control_points[3].y < 2.55, "Incorrect y coordinate of control point at position 3."
    assert 1.90 < spl.control_points[4].x < 1.91, "Incorrect x coordinate of control point at position 4."
    assert 2.98 < spl.control_points[4].y < 2.99, "Incorrect y coordinate of control point at position 4."
    assert -2.02 < spl.control_points[5].x < -2.01, "Incorrect x coordinate of control point at position 5."
    assert 3.52 < spl.control_points[5].y < 3.53, "Incorrect y coordinate of control point at position 5."
    assert 7.36 < spl.control_points[6].x < 7.37, "Incorrect x coordinate of control point at position 6."
    assert 0.90 < spl.control_points[6].y < 0.91, "Incorrect y coordinate of control point at position 6."
    assert 2.56 < spl.control_points[7].x < 2.57, "Incorrect x coordinate of control point at position 7."
    assert -1.17 < spl.control_points[7].y < -1.16, "Incorrect y coordinate of control point at position 7."
