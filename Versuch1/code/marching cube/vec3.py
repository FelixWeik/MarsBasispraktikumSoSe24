#!/usr/bin/python
import math


class Vec3:
    def __init__(self, a, b, c):
        self.x = a
        self.y = b
        self.z = c

    def __str__(self):
        return "(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")"

    def __add__(self, other):
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        return Vec3(self.x * other, self.y * other, self.z * other)

    # added this so we can calculate q
    def __truediv__(self, other):
        if other == 0:
            raise ValueError("Division by zero is not allowed")
        return Vec3(self.x / other, self.y / other, self.z / other)

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __ne__(self, other):
        return not self == other

    def __neg__(self):
        return -1 * self

    def __pos__(self):
        return self

    def dot(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z

    def __abs__(self):
        return math.sqrt(self.dot(self))

    def __hash__(self):
        return hash((self.x, self.y, self.z))


def vec_from_list(l):
    return Vec3(l[0], l[1], l[2])


if __name__ == "__main__":
    v = Vec3(4, 3, 2)
    w = Vec3(3, 1, .5)
    n = 3
    print(v)
    print(w)
    print(v + w)
    print(v - w)
    print(v * n)
    print(n * v)
    print(v == w)
    print(v != w)
    print(v.dot(w))
    print(abs(v))
