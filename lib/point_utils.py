from lib.vector import Vector
from shapely.geometry import Point

def add(point: Point, vector: Vector) -> Point:
    return Point(point.x + vector.x, point.y + vector.y)


def between_two_points(p1: Point, p2: Point, how_far: float) -> Point:
    if how_far > 1 or how_far < 0:
        raise Exception("factor between two points must be between 0 and 1")

    return Point(p1.x + how_far * (p2.x - p1.x), p1.y + how_far * (p2.y - p1.y))
