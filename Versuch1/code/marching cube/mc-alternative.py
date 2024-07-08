#!/usr/bin/python
from vec3 import Vec3, vec_from_list
import cube
from marchingviewer3d import Viewer3d
import multiprocessing


class Marching:
    # Add funtions to calculate the vertices and faces of the 3D object.
    # A 3D model consists of multiple polygons.  A polygon consists of vertices 
    # which are vec3 points and faces which are the indexes of the vertices.
    # A face of a polygon starts with the amount of points per polygon, for this task it will always be 3.
    #
    # For example a square made of two polygons could be:
    # vertices = [vec3(0, 0, 0), vec3(1, 0, 0), vec3(1, 1, 0), vec3(0, 1, 0)]
    # faces = [3, 0, 1, 2, 3, 3, 1, 2]
    # 4 vertices as corners for the square and the faces start with a three 
    # followed by three indexes corresponding to the vertice list
    # https://docs.pyvista.org/version/stable/examples/00-load/create-poly.html

    def __init__(self):
        self.vertices = cube.CubeVertices
        self.faces = []


    def get_CubeEdgeFlags(self):
        return cube.CubeEdgeFlags


    def get_CubeTriangles_from_cube_index(self, cube_index):
        return cube.CubeTriangles[cube_index]


    def get_cube_vertices_from_CubeEdges(self, edge_index):
        return cube.CubeEdges[edge_index]


    def get_CubeVertexDirections(self):
        return cube.CubeVertexDirections


    def get_CubeVertices_from_vertex_index(self, vertex_index):
        return cube.CubeVertices[vertex_index]


    def get_vertices_coord_from_triangles(self, cube_index, debug):
        # get the list of triangles of the cube
        cube_triangles = self.get_CubeTriangles_from_cube_index(cube_index)

        vertices_coord = []

        # for each value of the triangle list (values go 3 by 3)
        for i in range(0, len(cube_triangles), 3):
            # when we get the -1 value, it means that there are no triangles left
            if cube_triangles[i] == -1:
                # so we stop
                break

            # e1, e2, e3 represent the intersected edges
            # each group of three intersected edges form one triangle
            e1, e2, e3 = cube_triangles[i:i+3]
            if debug:
                print(f'e1: {e1}, e2: {e2}, e3: {e3}\n')

            # for each edge of the triangle
            for edge in (e1, e2, e3):
                if debug:
                    print(f'edge: {edge}')
                # get both vertices index that form the edge
                vertex_1_index, vertex_2_index = self.get_cube_vertices_from_CubeEdges(edge)
                if debug:
                    print(f'vertex_1_index: {vertex_1_index}')
                    print(f'vertex_2_index: {vertex_2_index}')

                # get coordinates of the vertices from their indexes
                vertex_1_coord = self.get_CubeVertices_from_vertex_index(vertex_1_index)
                vertex_2_coord = self.get_CubeVertices_from_vertex_index(vertex_2_index)
                vertices_coord.append([vertex_1_coord, vertex_2_coord])
                if debug:
                    print(f'vertex_1_coord: {vertex_1_coord}')
                    print(f'vertex_2_coord: {vertex_2_coord}\n')
            
        return vertices_coord


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
def march_cubes(data, debug):
    marching = data[0]
    length = data[1]
    subdiv = data[2]
    isoval = data[3]
    function = data[4]

    if debug:
        print('cube_index: 1\n')

    # get the coordinates of each vertex near the triangle's points
    vertices_coord = marching.get_vertices_coord_from_triangles(cube_index=1, debug=debug)

    if debug:
        print('\n')

    # for each pair of vertices that form the edge that is intersected by the triangle
    for (p1, p2) in vertices_coord:
        if debug:
            print(f'p1: {p1}, p2: {p2}')

        # transform the list of coordinates into a Vec3
        p1 = vec_from_list(p1)
        p2 = vec_from_list(p2)

        v1 = function(p1) - isoval
        v2 = function(p2) - isoval 
        if debug:
            print(f'v1: {v1}, v2: {v2}')

        # calculate q: a corner point of the triangle
        q = (v1*p2 - v2*p1) / (v1 - v2)
        if debug:
            print(f'q: {q}\n')

    return data


if __name__ == "__main__":
    if SHOW_SINGLE:
        mc = Marching()
        march_cubes((mc, LENGTH, SUBDIV, ISOVAL, sphere), debug=True)

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
