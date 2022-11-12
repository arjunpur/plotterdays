import random
from typing import Sequence 
from shapely.geometry import Point, Polygon

def generate_random_points_in_polygon(number: int, polygon: Polygon) -> Sequence[Point]:
    """
    Given a shapely.Polygon, generate a set of N random points in that polygon
    """
    points = []
    if len(polygon.bounds) == 0:
        raise Exception("polgon isn't bounded")
    minx, miny, maxx, maxy = polygon.bounds
    while len(points) < number:
        pnt = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
        if polygon.contains(pnt):
            points.append(pnt)
    return points
