from typing import Sequence
from lib.coordinates import pol2cart
from shapely.geometry import Point


class Spiral:
    def __init__(
        self,
        center: Point,
        start_radius: float,
        end_radius: float,
        space_per_loop: float,
        start_theta: float,
        theta_step: float,
    ):
        self.center = center
        self.start_radius = start_radius
        self.end_radius = end_radius
        self.space_per_loop = space_per_loop
        self.start_theta = start_theta
        self.theta_step = theta_step

    def get_points(self) -> Sequence[Point]:
        # This uses the basic formula for determining the radius to use
        # when drawing a point on a spiral (in polar coordinates)
        # The theta parameter will keep changing. `b` represents the space_per_loop
        # between each loop, and `a` represents where you're starting
        #
        # https://en.wikipedia.org/wiki/Archimedean_spiral
        def calculate_x_y(r: float, temp_theta: float) -> Point:
            point = pol2cart(r, temp_theta)
            centered = Point(point.x + self.center.x, point.y + self.center.y)
            return centered

        a = self.start_radius
        b = self.space_per_loop
        temp_theta = self.start_theta
        temp_radius = self.start_radius

        # Radius is updated by taking the starting radius of the circle and then
        # increasing it by a factor of `space_per_loop` scaled by the theta.
        # Scaling by the theta (which is incremented gradually) represents moving a point
        # in a circle at slowly increasing distances in polar coordinates
        r = a + (b * temp_theta)
        start_point = calculate_x_y(r, temp_theta)
        spiral_points = [start_point]

        while temp_radius < self.end_radius:
            temp_theta += self.theta_step
            r = a + (b * temp_theta)
            point = calculate_x_y(r, temp_theta)  # Calculate new radius after theta change
            temp_radius = r
            spiral_points.append(point)
        return spiral_points
