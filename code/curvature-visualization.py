#!/usr/bin/python

from cagd.spline import SplineSurface, Knots
from cagd.bezier import BezierSurface, BezierPatches
from cagd.vec import Vec3

f = lambda x, y: x * x / 5 + x * y / 4 + 10 / (1 + x * x + y * y) + y / 2
ctrl_pts = [[Vec3(x, y, f(x, y)) for x in range(-5, 5)] for y in range(-5, 5)]

d = 3
m = d + len(ctrl_pts) + 1
ku = Knots(m)
kv = Knots(m)
for i in range(d):
    ku[i] = 0
    kv[i] = 0
for i in range(m - d - d):
    ku[i + d] = i
    kv[i + d] = i
for i in range(d):
    ku[i + m - d] = m - d - d - 1
    kv[i + m - d] = m - d - d - 1

surface = SplineSurface((d, d))
surface.control_points = ctrl_pts
surface.knots = (ku, kv)

bezier_patches = surface.to_bezier_patches()
bezier_patches.refine(1)  # refine surface into more patches for more detailed coloring
bezier_patches.visualize_curvature(bezier_patches.CURVATURE_GAUSSIAN, bezier_patches.COLOR_MAP_LINEAR)

path = "surfaces.off"
f = open(path, 'w')
f.write(bezier_patches.export_standard_off())
f.close()
