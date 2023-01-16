import bezier
from lib.point_utils import add, between_two_points
from lib.shapely import generate_random_points_in_polygon
from lib.vector import Vector
import vsketch
import numpy as np
from typing import Tuple 
from shapely.geometry import LineString, Point, Polygon
from shapely.affinity import scale, translate


def create_curve_with_light_bend_and_noise(
        vsk: vsketch.Vsketch,
        start_end: Tuple[Point, Point],
        control_point_locations: Tuple[float, float],
        trunk_width: float,
        tip_width: float,
        control_point_precision: float,
        bend_clockwise: bool = True,
    ) -> np.ndarray:

    num_points = 1000

    start, end = start_end[0], start_end[1]
    noise_vals = vsk.noise([0., 1.], np.linspace(start.x, end.x, num_points) / 2)

    curve = create_curve_with_light_bend(
        start_end,
        control_point_locations,
        trunk_width,
        tip_width,
        control_point_precision,
        bend_clockwise,
    )
    
    points = curve.evaluate_multi(np.linspace(0, 1, num_points))
    new_points = points + noise_vals

    # Turn all tuples in zipped into a Point
    zipped = list(zip(new_points[0], new_points[1]))
    new_points  = [Point(p) for p in zipped]

    # Scale and translate the points to ensure the curve starts and ends
    # at (start, end)
    scale_factor_x = (end.x - start.x) / (new_points[-1].x - new_points[0].x)
    scale_factor_y = (end.y - start.y) / (new_points[-1].y - new_points[0].y)
    # translation_vector = Vector(sum(noise_vals[0]), sum(noise_vals[1]))

    start_point_displacement = Vector.from_two_points(new_points[0], start)
    end_point_displacement = Vector.from_two_points(new_points[-1], end)

    displacement_vector = start_point_displacement

    # zip the points in new_points together
    line = LineString(zipped)
    line = scale(line, xfact = scale_factor_x, yfact = scale_factor_y, origin = "centroid")
    line = translate(line, xoff = displacement_vector.x, yoff = displacement_vector.y)

    # TODO: Figure out why the line is still displaced by a bit 
    coords = line.coords
    x_coords = [coord[0] for coord in coords]
    y_coords = [coord[1] for coord in coords]

    return np.asarray([x_coords, y_coords])

def create_curve_with_light_bend(
        start_end: Tuple[Point, Point],
        control_point_locations: Tuple[float, float],
        trunk_width: float,
        tip_width: float,
        control_point_precision: float,
        bend_clockwise: bool = True,
    ) -> bezier.Curve:
    """
    Creates a curve that bends slightly.

    start_end: Where the curve starts and where it ends
    trajectory: What direction the curve is oriented towards
    control_point_locations: The points along the `trajectory` where we branch out to drop control points
    trunk_width: How far we branch out for the first control point
    tip_width: How far we branch out for the second control point
    control_point_precision: How random / precise our control point selection is - control points are chosen at random from within 
    bend_clockwise: Decides which direction to bend the curve. Clockwise `True` will use an 
    imaginary clock from the origin to determine the direction of the perpendicular vector
    a circle

    """
    start, end = start_end
    trajectory = Vector.from_two_points(start, end).normalize()

    first_point = between_two_points(start, end, control_point_locations[0]) 
    perpendicular = trajectory.get_perpendicular(clockwise = bend_clockwise).normalize()
    control_point_1_boundary = add(first_point, perpendicular * trunk_width)
    control_point_1_buffer = control_point_1_boundary.buffer(control_point_precision)

    second_point = between_two_points(start, end, control_point_locations[1]) 
    perpendicular = trajectory.get_perpendicular(clockwise = bend_clockwise).normalize()
    control_point_2_boundary = add(second_point, perpendicular * tip_width)
    control_point_2_buffer = control_point_2_boundary.buffer(control_point_precision)

    control_point_1 = generate_random_points_in_polygon(1, Polygon(control_point_1_buffer))[0] 
    control_point_2 = generate_random_points_in_polygon(1, Polygon(control_point_2_buffer))[0] 

    nodes = np.asarray([
        [start.x, control_point_1.x, control_point_2.x, end.x],
        [start.y, control_point_1.y, control_point_2.y, end.y]
    ])
    return bezier.Curve(
        nodes,
        degree = 3
    )

def draw_cubic_bezier(vsk: vsketch.Vsketch, curve: bezier.Curve):
    nodes = curve.nodes
    assert(nodes.shape == (2, 4)) # We expect a cubic Bezier curve
    vsk.bezier(
        nodes[0][0],
        nodes[1][0],
        nodes[0][1],
        nodes[1][1],
        nodes[0][2],
        nodes[1][2],
        nodes[0][3],
        nodes[1][3] 
    )
