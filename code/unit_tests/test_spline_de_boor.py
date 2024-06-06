from cagd.spline import Knots, Spline
from cagd.vec import Vec2


def spline_config1():
    kts = Knots(9)
    kts.knots = [0, 1, 2, 3, 4, 5, 6, 7]
    degree = 3
    control_points = [Vec2(0, 0), Vec2(0.5, 4), Vec2(4, 4.5), Vec2(5, 1)]
    spl = Spline(degree)
    spl.knots = kts
    spl.control_points = control_points
    return spl


def spline_config2():
    kts = Knots(11)
    kts.knots = [0, 0, 2.5, 3, 4.5, 5, 6.2, 7, 7]
    degree = 3
    control_points = [Vec2(0, 0), Vec2(0.5, 4), Vec2(4, 4.5), Vec2(5, 1), Vec2(4, 8)]
    spl = Spline(degree)
    spl.knots = kts
    spl.control_points = control_points
    return spl


def test_correct_return_type():
    spl = spline_config1()

    p = spl.de_boor(3, 1)
    assert type(p) == list, "De Boor should return a list."
    assert type(p[0]) == Vec2, "The returned list should contain cagd.Vec2 type elements."


def test_only_last_point_config1():
    spl = spline_config1()

    points = spl.de_boor(3, 1)
    assert len(points) == 1, "Wrong number of points returned."
    p = points[0]
    assert 0.99 < p.x < 1.0, "x coordinate is incorrect."
    assert 3.41 < p.y < 3.42, "y coordinate is incorrect."


def test_just_return_control_points_config1():
    spl = spline_config1()

    points = spl.de_boor(3, 4)
    assert len(points) == 4, "Wrong number of points returned."
    for c, p in zip(spl.control_points, points):
        assert c == p, "Returned points are not the spline's control points."


def test_return_first_step_config1():
    spl = spline_config1()

    points = spl.de_boor(3, 3)
    assert len(points) == 3, "Wrong number of points returned."
    assert 0.33 < points[0].x < 0.34, "Point d10 has wrong x coordinate."
    assert 2.66 < points[0].y < 2.67, "Point d10 has wrong y coordinate."
    assert 1.66 < points[1].x < 1.67, "Point d11 has wrong x coordinate."
    assert 4.16 < points[1].y < 4.17, "Point d11 has wrong y coordinate."
    assert 3.99 < points[2].x < 4.01, "Point d12 has wrong x coordinate."
    assert 4.49 < points[2].y < 4.51, "Point d12 has wrong y coordinate."


def test_return_second_step_config1():
    spl = spline_config1()

    points = spl.de_boor(3, 2)
    assert len(points) == 2, "Wrong number of points returned."
    assert 0.99 < points[0].x < 1.0, "Point d20 has wrong x coordinate."
    assert 3.41 < points[0].y < 3.42, "Point d20 has wrong y coordinate."
    assert 1.66 < points[1].x < 1.67, "Point d21 has wrong x coordinate."
    assert 4.16 < points[1].y < 4.17, "Point d21 has wrong y coordinate."


def test_t_outside_of_support_config1_upper():
    spl = spline_config1()
    try:
        spl.de_boor(5, 1)
    except AssertionError:
        assert True
        return
    assert False, "t=5 is not in the support and de Boor should not be executed."


def test_t_outside_of_support_lower():
    spl = spline_config1()
    try:
        spl.de_boor(2, 1)
    except AssertionError:
        assert True
        return
    assert False, "t=5 is not in the support and de Boor should not be executed."


def test_t_outside_of_support_config2_upper():
    spl = spline_config2()
    try:
        spl.de_boor(6, 1)
    except AssertionError:
        assert True
        return
    assert False, "t=5 is not in the support and de Boor should not be executed."


def test_t_outside_of_support_config2_lower():
    spl = spline_config2()
    try:
        spl.de_boor(2, 1)
    except AssertionError:
        assert True
        return
    assert False, "t=5 is not in the support and de Boor should not be executed."


def test_only_last_point_config2():
    spl = spline_config2()

    points = spl.de_boor(3, 1)
    assert len(points) == 1, "Wrong number of points returned."
    p = points[0]
    print(p.x, p.y)
    assert 0.54 < p.x < 0.56, "x coordinate is incorrect."
    assert 3.02 < p.y < 3.03, "y coordinate is incorrect."


def test_just_return_control_points_config2():
    spl = spline_config2()

    points = spl.de_boor(3, 4)
    assert len(points) == 4, "Wrong number of points returned."
    for c, p in zip(spl.control_points, points):
        assert c == p, "Returned points are not the spline's control points."


def test_return_first_step_config2():
    spl = spline_config2()

    points = spl.de_boor(3, 3)
    assert len(points) == 3, "Wrong number of points returned."
    assert 0.33 < points[0].x < 0.34, "Point d10 has wrong x coordinate."
    assert 2.66 < points[0].y < 2.67, "Point d10 has wrong y coordinate."
    assert 1.19 < points[1].x < 1.21, "Point d11 has wrong x coordinate."
    assert 4.09 < points[1].y < 4.11, "Point d11 has wrong y coordinate."
    assert 3.99 < points[2].x < 4.01, "Point d12 has wrong x coordinate."
    assert 4.49 < points[2].y < 4.51, "Point d12 has wrong y coordinate."


def test_return_second_step_config2():
    spl = spline_config2()

    points = spl.de_boor(3, 2)
    assert len(points) == 2, "Wrong number of points returned."
    assert 0.54 < points[0].x < 0.56, "Point d20 has wrong x coordinate."
    assert 3.02 < points[0].y < 3.03, "Point d20 has wrong y coordinate."
    assert 1.19 < points[1].x < 1.21, "Point d21 has wrong x coordinate."
    assert 4.09 < points[1].y < 4.11, "Point d21 has wrong y coordinate."


def test_only_last_point_config2_second_segment():
    spl = spline_config2()

    points = spl.de_boor(4, 1)
    assert len(points) == 1, "Wrong number of points returned."
    p = points[0]
    assert 2.99 < p.x < 3.0, "x coordinate is incorrect."
    assert 3.94 < p.y < 3.95, "y coordinate is incorrect."


def test_just_return_control_points_config2_second_segment():
    spl = spline_config2()

    points = spl.de_boor(4, 4)
    assert len(points) == 4, "Wrong number of points returned."
    for c, p in zip(spl.control_points, points):
        assert c == p, "Returned points are not the spline's control points."


def test_return_first_step_config2_second_segment():
    spl = spline_config2()

    points = spl.de_boor(4, 3)
    assert len(points) == 3, "Wrong number of points returned."
    assert 0.44 < points[0].x < 0.45, "Point d10 has wrong x coordinate."
    assert 3.55 < points[0].y < 3.56, "Point d10 has wrong y coordinate."
    assert 2.59 < points[1].x < 2.61, "Point d11 has wrong x coordinate."
    assert 4.29 < points[1].y < 4.31, "Point d11 has wrong y coordinate."
    assert 4.31 < points[2].x < 4.32, "Point d12 has wrong x coordinate."
    assert 3.40 < points[2].y < 3.41, "Point d12 has wrong y coordinate."


def test_return_second_step_config2_second_segment():
    spl = spline_config2()

    points = spl.de_boor(4, 2)
    assert len(points) == 2, "Wrong number of points returned."
    assert 2.06 < points[0].x < 2.07, "Point d20 has wrong x coordinate."
    assert 4.11 < points[0].y < 4.12, "Point d20 has wrong y coordinate."
    assert 3.45 < points[1].x < 3.46, "Point d21 has wrong x coordinate."
    assert 3.85 < points[1].y < 3.86, "Point d21 has wrong y coordinate."
