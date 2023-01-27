from lib.random_utils import random_n_elements_across_k_partitions
import vsketch
import numpy as np
from random import choice
from lib.curves import create_curve_with_light_bend, create_curve_with_light_bend_and_noise, draw_cubic_bezier
from lib.grid import Grid, create_grid_with_padding, draw_grid
from lib.point_utils import add
from lib.vector import Vector, draw_vector, rotate_vector
from shapely.geometry import Point

class RecursiveVeinsSketch(vsketch.SketchClass):

    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("9in", "11in", landscape=True, center=False)
        vsk.scale("in")

        vsk.stroke(1)
        grid = create_grid_with_padding(padding = 0.5, cell_size = (0.50, 0.50)) 
        rect_cell = grid.get_cell_at_index(0, 0)
        sub_grid = Grid(rect_cell.width * 8, rect_cell.height * 8, rect_cell.top_left, 10, 10) 
        end_cell = sub_grid.get_cell_at_index(0, 8)
        start_cell = sub_grid.get_cell_at_index(5, 9)
        
        vsk.point(start_cell.center.x, start_cell.center.y)
        vsk.point(end_cell.center.x, end_cell.center.y)

        draw_grid(vsk, sub_grid)

        unit_trajectory = Vector.from_two_points(start_cell.center, end_cell.center).normalize()
        length = Vector.from_two_points(start_cell.center, end_cell.center).length / 2

        vsk.stroke(2)
        self.create_branched_curves(
            vsk,
            start_cell.center,
            length, 
            unit_trajectory,
            count = 0,
            trunk_width = 0.2,
            control_point_precision = 0.01,
        ) 

    def create_branched_curves(
        self,
        vsk: vsketch.Vsketch,
        start_point: Point,
        length: float,
        unit_trajectory: Vector,
        count: int,
        trunk_width: float,
        control_point_precision: float,
    ):
        if count > 2:
            return

        # Draw the curve as described by the given parameters 
        bend_clockwise = unit_trajectory.x > 0
        end_point = add(start_point, unit_trajectory * length)

        # With noise
        # curve = create_curve_with_light_bend_and_noise(
        #     vsk,
        #     start_end = (start_point, end_point),
        #     control_point_locations = (0.25, 0.75),
        #     trunk_width = trunk_width,
        #     tip_width = 0.0,
        #     control_point_precision = 0.01,
        #     bend_clockwise=bend_clockwise,
        # )
        # vsk.polygon(curve[0], curve[1])
        # indices = random_n_elements_across_k_partitions(curve[0], 3, 3)
        # start_points = [curve[0][indices], curve[1][indices]]
        
        # No Noise approach
        no_noise_curve = create_curve_with_light_bend(
            start_end = (start_point, end_point),
            control_point_locations = (0.25, 0.75),
            trunk_width = trunk_width,
            tip_width = 0.0,
            control_point_precision = 0.01,
            bend_clockwise=bend_clockwise,
        )
        draw_cubic_bezier(vsk, no_noise_curve)
        points = no_noise_curve.evaluate_multi(np.linspace(0, 1, 1000))
        indices = random_n_elements_across_k_partitions(points[0], 3, 3)
        start_points = [points[0][indices], points[1][indices]]


        start_angle = choice([30])
        for i in range(len(start_points[0])):
            # Pick a random angle to rotate
            angle = start_angle - (i * 5) 
            new_length = length * 0.5
            new_count = count + 1
            rotated_trajectory = rotate_vector(unit_trajectory, angle, degrees = True)
            self.create_branched_curves(
                vsk,
                Point(start_points[0][i], start_points[1][i]),
                new_length,
                rotated_trajectory,
                new_count,
                trunk_width * 0.6,
                control_point_precision * 0.6
            )
