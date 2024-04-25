#!/usr/bin/python
import math


class Vec3:
    def __init__(self, a, b, c):
        self.x = a
        self.y = b
        self.z = c

    def __str__(self):
        return "(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")"

    def __repr__(self):
        return self.__str__()

    def __unicode__(self):
        return 'u' + self.__str__()

    def __add__(self, other):
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        return Vec3(self.x * other, self.y * other, self.z * other)

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

    def __copy__(self):
        return Vec3(self.x, self.y, self.z)


class Vec2:
    def __init__(self, a, b):
        self.x = a
        self.y = b

    def __str__(self):
        return "(" + str(self.x) + "," + str(self.y) + ")"

    def __repr__(self):
        return self.__str__()

    def __unicode__(self):
        return 'u' + self.__str__()

    def __add__(self, other):
        return Vec2(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Vec2(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        return Vec2(self.x * other, self.y * other)

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return not self == other

    def __neg__(self):
        return -1 * self

    def __pos__(self):
        return self

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def __abs__(self):
        return math.sqrt(self.dot(self))

    def __hash__(self):
        return hash((self.x, self.y))

    def __copy__(self):
        return Vec2(self.x, self.y)
