#!/usr/bin/python

from cagd.polyline import Polyline
from cagd.spline import Spline, SplineSurface, Knots
from cagd.bezier import BezierSurface, BezierPatches
from cagd.vec import Vec2, Vec3
from cagd.viewer3d import Viewer3d

import cagd.scene_2d as scene_2d

pts = [Vec2(0.05, 5.5),
       Vec2(1.5, 5),
       Vec2(2, 4),
       Vec2(1.7, 2.5),
       Vec2(0.7, 1.8),
       Vec2(2, 1.3),
       Vec2(2, 0.9),
       Vec2(1.2, 0.8),
       Vec2(.7, 0.4),
       Vec2(.7, -1),
       Vec2(.7, -2.8),
       Vec2(2, -4),
       Vec2(2, -4.6), ]

spl = Spline.interpolate_cubic(Spline.INTERPOLATION_CHORDAL, pts, Knots(1))
spl.set_color("#0000ff")
sc = scene_2d.Scene()
sc.set_resolution(900)
sc.add_element(spl)


surface = spl.generate_rotation_surface(8)

print(surface.knots)

v = Viewer3d()
# show control points of rotated surface
cps = [pt for pts in surface.control_points for pt in pts]
v.display_points(cps, Vec3(0, 0, 0), "red")

bezier_patches = surface.to_bezier_patches()
# show points of bezier patches
cps = [pt for pts in surface.control_points for pt in pts]
v.display_points(cps, Vec3(0,0,0), "green")
v.display_object(bezier_patches, Vec3(0,0,0))

bezier_patches.refine(2)
# show refined object
v.display_object(bezier_patches, Vec3(-5,5,0))
v.show()
