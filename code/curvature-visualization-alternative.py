#!/usr/bin/python

from cagd.spline import Spline, SplineSurface, Knots
from cagd.bezier import BezierSurface, BezierPatches
from cagd.vec import Vec2, Vec3
from cagd.viewer3d import Viewer3d
import multiprocessing
import os

# if True only show a single object otherwise show all objects
show_single = True
# if True use figure from task 3 otherwise use deformed plane
figure = False


# parallel compute curvature visualization
def visualize_curvature(patches):
    patches[0].visualize_curvature(patches[1], patches[2])
    return patches


def get_surface():
    if figure:
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
        rot_surface = spl.generate_rotation_surface(6)
    else:
        f = lambda x, y: x * x / 5 + x * y / 4 + 10 / (1 + x * x + y * y) + y / 2
        ctrl_pts = [[Vec3(x, y, f(x, y)) for x in range(-5, 6)] for y in range(5, -6, -1)]

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

        rot_surface = SplineSurface((d, d))
        rot_surface.control_points = ctrl_pts
        rot_surface.knots = (ku, kv)

    return rot_surface


if __name__ == "__main__":
    surface = get_surface()

    v = Viewer3d()
    bezier_patches = surface.to_bezier_patches()
    bezier_patches.refine(1)

    if not show_single:
        v.add_text("                                        CURVATURE_GAUSSIAN          CURVATURE_AVERAGE          CURVATURE_PRINCIPAL_MAX          CURVATURE_PRINCIPAL_MIN \
            \n\n\n\n\n\n\n\nCOLOR_MAP_LINEAR\n\n\n\n\n\n\n\nCOLOR_MAP_CUT\n\n\n\n\n\n\n\nCOLOR_MAP_CLASSIFICATION")

        threads = []
        bpatches = []
        indexes = []
        inputData = []
        for i in range(4):
            for j in range(4, 7):
                inputData.append((bezier_patches, i, j))

        p = multiprocessing.Pool(os.cpu_count())
        data = p.map(visualize_curvature, inputData)

        for di in data:
            i = di[1]
            j = di[2]
            v.display_object(di[0], Vec3(-10 * i + 10 * j, 10 * i + 10 * j, 0))

        p.close()
        v.show()

    else:
        bezier_patches.visualize_curvature(bezier_patches.CURVATURE_GAUSSIAN, bezier_patches.COLOR_MAP_LINEAR)

        v.display_object(bezier_patches, Vec3(0, 0, 0))
        v.show()
