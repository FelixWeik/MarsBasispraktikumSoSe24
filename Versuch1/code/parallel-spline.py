#!/usr/bin/python

from cagd.polyline import Polyline
from cagd.spline import Spline, Knots
from cagd.vec import Vec2
import cagd.scene_2d as scene_2d

pts = [Vec2(0, .4), Vec2(.8, .8), Vec2(.5, 1.2), Vec2(-.03, .4), Vec2(.4, 0), Vec2(1, .2)]
s1 = Spline.interpolate_cubic(Spline.INTERPOLATION_CHORDAL, pts, Knots(1))

with_added = Spline.interpolate_cubic(Spline.INTERPOLATION_CHORDAL, pts, Knots(1))

pts = [Vec2(0, 2.5), Vec2(-1, 1), Vec2(1, -1), Vec2(0, -2.5), Vec2(-1, -1), Vec2(1, 1)]
s1 = Spline.interpolate_cubic_periodic(pts)

s1.set_color("#0000ff")

sc = scene_2d.Scene()
sc.set_resolution(900)


with_added.insert_knot(0.2)
with_added.insert_knot(1.2)
with_added.insert_knot(1.8)
with_added.insert_knot(2.2)
with_added.insert_knot(3.2)
with_added.set_color("#ff0000")

print(with_added.tangent(0.4))

sc.add_element(with_added)


sc.add_element(s1)

for i in [-1, 1]:
    para = s1.generate_parallel(i * 0.025, 0.005)
    para.set_color("#999999")
    sc.add_element(para)

sc.write_image()
sc.show()
