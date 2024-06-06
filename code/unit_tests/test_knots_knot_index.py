from cagd.spline import Knots


def config1():
    kts = Knots(6)
    kts.knots = [0, 2, 3, 5, 7, 8]
    return kts


def config2():
    kts = Knots(8)
    kts.knots = [0, 0, 1, 1, 1, 2, 2]
    return kts


def test_regular_intervals():
    kts = config1()

    assert kts.knot_index(-0.1) is None
    assert kts.knot_index(0) == 0
    assert kts.knot_index(0.2) == 0
    assert kts.knot_index(1) == 0
    assert kts.knot_index(2) == 1
    assert kts.knot_index(4) == 2
    assert kts.knot_index(7) == 4
    """
    Theoretically, t=8 is not inside the spline anymore,
    but we allow it so that the spline is easier to draw.
    """
    assert kts.knot_index(8) == 4, "Edge case, see comment in code."
    assert kts.knot_index(8.1) is None


def test_multi_knots():
    kts = config2()

    assert kts.knot_index(-0.1) is None
    assert kts.knot_index(0) == 1
    assert kts.knot_index(0.1) == 1
    assert kts.knot_index(1) == 4
    assert kts.knot_index(1.1) == 4
    assert kts.knot_index(2) == 4
    assert kts.knot_index(2.1) is None
