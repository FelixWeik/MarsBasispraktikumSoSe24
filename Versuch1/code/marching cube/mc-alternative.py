#!/usr/bin/python
from vec3 import Vec3, vec_from_list
import cube
from marchingviewer3d import Viewer3d
import multiprocessing


class Marching:
    # Add funtions to calculate the vertices and faces of the 3D object.
    # A 3D model consists of multiple polygons.  A polygon consists of vertices 
    # which are vec3 points and faces which are the indexes of the vertices.
    # A face of a polygon starts with the amount of points per polygon, for this task it wil always be 3.
    #
    # For example a square made of two polygons could be:
    # vertices = [vec3(0, 0, 0), vec3(1, 0, 0), vec3(1, 1, 0), vec3(0, 1, 0)]
    # faces = [3, 0, 1, 2, 3, 3, 1, 2]
    # 4 vertices as corners for the square and the faces start with a three 
    # followed by three indexes corresponding to the vertice list
    # https://docs.pyvista.org/version/stable/examples/00-load/create-poly.html

    def __init__(self):
        self.vertices = []
        self.faces = []


LENGTH = 1
ISOVAL = 1  # iso-value
SUBDIV = 32  # resolution
SHOW_SINGLE = True


def sphere(vec):
    return abs(vec)


def octahedron(vec):
    return abs(vec.x) + abs(vec.y) + abs(vec.z)


def cube_func(vec):
    return max(max(abs(vec.x), abs(vec.y)), abs(vec.z))


def torus(vec):
    x = vec.x
    y = vec.y
    z = vec.z
    c = x * x + y * y + z * z + .7 * .7 - .2 * .2
    d = 4 * .7 * .7 * (x * x + y * y)
    return c * c - d


# calculate vertices and faces of marching
def march_cubes(data):
    marching = data[0]
    length = data[1]
    subdiv = data[2]
    isoval = data[3]
    function = data[4]

    cube_size = (length * 2) / subdiv

    for x in range(0, SUBDIV):
        for y in range(0, SUBDIV):
            for z in range(0, SUBDIV):
                cube_start = Vec3(-length + cube_size * x, -length + cube_size * y, -length + cube_size * z)

                corner_index = 0

                intersections = [-1] * 12

                for corner in range(0, 8):
                    rel_corner = Vec3(cube.CubeVertices[corner][0], cube.CubeVertices[corner][1], cube.CubeVertices[corner][2])
                    corner_pos = rel_corner * cube_size + cube_start
                    if function(corner_pos) - isoval < 0:
                        corner_index += 2 ** corner

                edge_index = cube.CubeEdgeFlags[corner_index]

                for edge in range(0, 12):
                    if edge_index & (2 ** edge) > 0:
                        [corner1, corner2] = cube.CubeEdges[edge]
                        rel_corner1 = Vec3(cube.CubeVertices[corner1][0], cube.CubeVertices[corner1][1],
                                           cube.CubeVertices[corner1][2])
                        corner_pos1 = rel_corner1 * cube_size + cube_start
                        value1 = function(corner_pos1) - isoval
                        rel_corner2 = Vec3(cube.CubeVertices[corner2][0], cube.CubeVertices[corner2][1],
                                           cube.CubeVertices[corner2][2])
                        corner_pos2 = rel_corner2 * cube_size + cube_start
                        value2 = function(corner_pos2) - isoval
                        intersections[edge] = (value1 * corner_pos2 - value2 * corner_pos1) * (1.0 / (value1 - value2))

                edge_points = cube.CubeTriangles[corner_index]

                for vertex in range(0, 16, 3):
                    if edge_points[vertex] >= 0:
                        marching.faces.append(3)

                        marching.faces.append(len(marching.vertices))
                        marching.vertices.append(intersections[edge_points[vertex]])

                        marching.faces.append(len(marching.vertices))
                        marching.vertices.append(intersections[edge_points[vertex + 1]])

                        marching.faces.append(len(marching.vertices))
                        marching.vertices.append(intersections[edge_points[vertex + 2]])

    return data


if __name__ == "__main__":
    if SHOW_SINGLE:
        mc = Marching()
        march_cubes((mc, LENGTH, SUBDIV, ISOVAL, sphere))

        v = Viewer3d()
        v.display_marching_cube(mc, Vec3(0, 0, 0))
        v.show()
    else:
        mc = Marching()
        inputData = [(mc, LENGTH, SUBDIV, ISOVAL, cube_func),
                     (mc, LENGTH, SUBDIV, ISOVAL, sphere),
                     (mc, LENGTH, SUBDIV, ISOVAL, octahedron),
                     (mc, LENGTH, SUBDIV, 0, torus)]

        p = multiprocessing.Pool(4)
        data = p.map(march_cubes, inputData)

        v = Viewer3d()
        v.display_marching_cube(data[0][0], Vec3(0, 0, 0))
        v.display_marching_cube(data[1][0], Vec3(-2, 2, 0))
        v.display_marching_cube(data[2][0], Vec3(2, 2, 0))
        v.display_marching_cube(data[3][0], Vec3(0, 4, 0))

        v.show()
