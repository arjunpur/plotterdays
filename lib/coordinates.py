import numpy as np
from shapely.geometry import Point


def cart2pol(x: float, y: float) -> Point:
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return Point(rho, phi)


def pol2cart(rho: float, phi: float) -> Point:
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return Point(x, y)
