#!/usr/bin/python

from cagd.polyline import Polyline
from cagd.spline import Spline, SplineSurface, Knots
from cagd.bezier import BezierSurface, BezierPatches
from cagd.vec import Vec2, Vec3
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
# you can activate these lines to view the input spline
# spl.set_color("#0000ff")
# sc = scene_2d.Scene()
# sc.set_resolution(900)
# sc.add_element(spl)
# sc.write_image()
# sc.show()

surface = spl.generate_rotation_surface(6)

bezier_patches = surface.to_bezier_patches()

bezier_patches.refine(2)
path = "surfaces.off"
f = open(path, 'w')
f.write(bezier_patches.export_standard_off())
f.close()
