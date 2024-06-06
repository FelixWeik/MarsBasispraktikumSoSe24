from cagd.spline import Spline
from cagd.vec import Vec2


def config1():
    pts = [Vec2(0.2, 3), Vec2(5, 1), Vec2(3, -0.2), Vec2(1, 2), Vec2(1, 3)]
    return pts


def test_interpolate_equidistant_degree():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_EQUIDISTANT, points)
    assert spl.degree == 3, "Incorrect degree."


def test_interpolate_equidistant_knot_amount():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_EQUIDISTANT, points)
    kts = [0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4]
    assert len(spl.knots) == len(kts), "Incorrect number of Knots."


def test_interpolate_equidistant_knot_values():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_EQUIDISTANT, points)
    kts = [0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4]
    assert len(spl.knots) == len(kts), "Incorrect number of Knots."
    for n, (i, j) in enumerate(zip(spl.knots, kts)):
        assert i == j, f"Wrong knot value at pos {n}."


def test_interpolate_equidistant_control_points_amount():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_EQUIDISTANT, points)
    assert len(spl.control_points) == 7, "Incorrect number of control points."


def test_interpolate_equidistant_control_points_values():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_EQUIDISTANT, points)
    assert 0.19 < spl.control_points[0].x < 0.21, "Incorrect x coordinate of control point at position 0."
    assert 2.99 < spl.control_points[0].y < 3.01, "Incorrect y coordinate of control point at position 0."
    assert 2.39 < spl.control_points[1].x < 2.4, "Incorrect x coordinate of control point at position 1."
    assert 2.34 < spl.control_points[1].y < 2.36, "Incorrect y coordinate of control point at position 1."
    assert 6.78 < spl.control_points[2].x < 6.79, "Incorrect x coordinate of control point at position 2."
    assert 1.04 < spl.control_points[2].y < 1.05, "Incorrect y coordinate of control point at position 2."
    assert 2.65 < spl.control_points[3].x < 2.66, "Incorrect x coordinate of control point at position 3."
    assert -1.20 < spl.control_points[3].y < -1.19, "Incorrect y coordinate of control point at position 3."
    assert 0.58 < spl.control_points[4].x < 0.59, "Incorrect x coordinate of control point at position 4."
    assert 2.54 < spl.control_points[4].y < 2.56, "Incorrect y coordinate of control point at position 4."
    assert 0.86 < spl.control_points[5].x < 0.87, "Incorrect x coordinate of control point at position 5."
    assert 2.84 < spl.control_points[5].y < 2.85, "Incorrect y coordinate of control point at position 5."
    assert 0.99 < spl.control_points[6].x < 1.01, "Incorrect x coordinate of control point at position 6."
    assert 2.99 < spl.control_points[6].y < 3.01, "Incorrect y coordinate of control point at position 6."


def test_interpolate_chordal_knot_amount():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_CHORDAL, points)
    assert len(spl.knots) == 11, "Incorrect number of Knots."


def test_interpolate_chordal_knot_values():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_CHORDAL, points)
    assert len(spl.knots) == 11, "Incorrect number of Knots."
    for i in range(4):
        assert spl.knots[i] == 0, f"Wrong knot value at pos {i}."
    assert 5.19 < spl.knots[4] < 5.21, "Wrong knot value at pos 5."
    assert 7.53 < spl.knots[5] < 7.54, "Wrong knot value at pos 6."
    assert 10.50 < spl.knots[6] < 10.51, "Wrong knot value at pos 7."
    for i in range(4):
        assert 11.50 < spl.knots[i + 7] < 11.51, f"Wrong knot value at pos {i + 7}."


def test_interpolate_chordal_control_points_amount():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_CHORDAL, points)
    assert len(spl.control_points) == 7, "Incorrect number of control points."


def test_interpolate_chordal_control_points_values():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_CHORDAL, points)
    assert 0.19 < spl.control_points[0].x < 0.21, "Incorrect x coordinate of control point at position 0."
    assert 2.99 < spl.control_points[0].y < 3.01, "Incorrect y coordinate of control point at position 0."
    assert 2.89 < spl.control_points[1].x < 2.9, "Incorrect x coordinate of control point at position 1."
    assert 2.59 < spl.control_points[1].y < 2.60, "Incorrect y coordinate of control point at position 1."
    assert 6.80 < spl.control_points[2].x < 6.81, "Incorrect x coordinate of control point at position 2."
    assert 1.99 < spl.control_points[2].y < 2.00, "Incorrect y coordinate of control point at position 2."
    assert 2.62 < spl.control_points[3].x < 2.63, "Incorrect x coordinate of control point at position 3."
    assert -1.09 < spl.control_points[3].y < -1.08, "Incorrect y coordinate of control point at position 3."
    assert 0.87 < spl.control_points[4].x < 0.88, "Incorrect x coordinate of control point at position 4."
    assert 1.36 < spl.control_points[4].y < 1.37, "Incorrect y coordinate of control point at position 4."
    assert 0.97 < spl.control_points[5].x < 0.98, "Incorrect x coordinate of control point at position 5."
    assert 2.67 < spl.control_points[5].y < 2.68, "Incorrect y coordinate of control point at position 5."
    assert 0.99 < spl.control_points[6].x < 1.01, "Incorrect x coordinate of control point at position 6."
    assert 2.99 < spl.control_points[6].y < 3.01, "Incorrect y coordinate of control point at position 6."


def test_interpolate_centripetal_knot_amount():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_CENTRIPETAL, points)
    assert len(spl.knots) == 11, "Incorrect number of Knots."


def test_interpolate_centripetal_knot_values():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_CENTRIPETAL, points)
    assert len(spl.knots) == 11, "Incorrect number of Knots."
    for i in range(4):
        assert spl.knots[i] == 0, f"Wrong knot value at pos {i}."
    assert 2.28 < spl.knots[4] < 2.29, "Wrong knot value at pos 5."
    assert 3.80 < spl.knots[5] < 3.81, "Wrong knot value at pos 6."
    assert 5.53 < spl.knots[6] < 5.54, "Wrong knot value at pos 7."
    for i in range(4):
        assert 6.53 < spl.knots[i + 7] < 6.54, f"Wrong knot value at pos {i + 7}."


def test_interpolate_centripetal_control_points_values():
    points = config1()
    spl = Spline.interpolate_cubic(Spline.INTERPOLATION_CENTRIPETAL, points)
    assert 0.19 < spl.control_points[0].x < 0.21, "Incorrect x coordinate of control point at position 0."
    assert 2.99 < spl.control_points[0].y < 3.01, "Incorrect y coordinate of control point at position 0."
    assert 2.6 < spl.control_points[1].x < 2.61, "Incorrect x coordinate of control point at position 1."
    assert 2.44 < spl.control_points[1].y < 2.45, "Incorrect y coordinate of control point at position 1."
    assert 6.62 < spl.control_points[2].x < 6.63, "Incorrect x coordinate of control point at position 2."
    assert 1.5 < spl.control_points[2].y < 1.52, "Incorrect y coordinate of control point at position 2."
    assert 2.66 < spl.control_points[3].x < 2.67, "Incorrect x coordinate of control point at position 3."
    assert -1.18 < spl.control_points[3].y < -1.17, "Incorrect y coordinate of control point at position 3."
    assert 0.76 < spl.control_points[4].x < 0.77, "Incorrect x coordinate of control point at position 4."
    assert 1.96 < spl.control_points[4].y < 1.97, "Incorrect y coordinate of control point at position 4."
    assert 0.93 < spl.control_points[5].x < 0.94, "Incorrect x coordinate of control point at position 5."
    assert 2.72 < spl.control_points[5].y < 2.73, "Incorrect y coordinate of control point at position 5."
    assert 0.99 < spl.control_points[6].x < 1.01, "Incorrect x coordinate of control point at position 6."
    assert 2.99 < spl.control_points[6].y < 3.01, "Incorrect y coordinate of control point at position 6."
