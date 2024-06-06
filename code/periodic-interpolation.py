#!/usr/bin/python
import math

import cagd.scene_2d as scene
from cagd.vec import Vec2
from cagd.spline import Spline, Knots
from cagd.polyline import Polyline
from test_spline_interpolate_cubic_periodic import *
from test_utils_solve_almost_tridiagonal_equation import test_solve_almost_tridiagonal_equation


# returns a list of num_samples points that are uniformly distributed on the unit circle
def unit_circle_points(num_samples):

    res = []
    for i in range(num_samples):
        rel = i / num_samples  # Calculate the normalized relative position
        angle = 2 * math.pi * rel  # Compute the angle in radians
        x = math.cos(angle)  # X-coordinate on the unit circle
        y = math.sin(angle)  # Y-coordinate on the unit circle
        res.append(Vec2(x, y))  # Append the point to the result list

    return res


# calculates the deviation between the given spline and a unit circle
def calculate_circle_deviation(spline):
    samples = 100
    accError = 0
    maxError = 0
    spline_knot_len = len(spline.knots.knots) - 1 - 2 * spline.degree
    for i in range(samples):
        rel = i / samples
        rel_circ = (i + 1) / samples
        circle = Vec2(math.cos(2 * math.pi * rel_circ), math.sin(2 * math.pi * rel_circ))
        circle_len = math.sqrt(circle.dot(circle))
        spline_vec = spline.evaluate(rel * spline_knot_len + spline.degree)
        spline_len = math.sqrt(spline_vec.dot(spline_vec))
        error = circle_len - spline_len
        if math.fabs(error) > math.fabs(maxError):
            maxError = error
        accError += error
    medianError = accError / samples

    print("Max Error: ", maxError)
    print("Median Error: ", medianError)




test_solve_almost_tridiagonal_equation()
test_interpolate_periodic_degree()
test_interpolate_periodic_knot_amount()
test_interpolate_periodic_knot_values()
test_interpolate_periodic_control_points_amount()
test_interpolate_periodic_control_points_values()

# interpolate 6 points with a periodic spline to create the number "8"
pts = [Vec2(0, 2.5), Vec2(-1, 1), Vec2(1, -1), Vec2(0, -2.5), Vec2(-1, -1), Vec2(1, 1)]
pts_line = Polyline()
pts_line.points = pts
pts_line.set_color("red")
s = Spline.interpolate_cubic_periodic(pts)
p = s.get_polyline_from_control_points()
p.set_color("blue")
sc = scene.Scene()
sc.set_resolution(900)
sc.add_element(s)
sc.add_element(p)

# generate a spline that approximates the unit circle
n = 8
circle_pts = unit_circle_points(n)
circle = Spline.interpolate_cubic_periodic(circle_pts)
sc.add_element(circle)
calculate_circle_deviation(circle)

sc.write_image()
sc.show()
