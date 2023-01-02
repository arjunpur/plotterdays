import vsketch
import numpy as np
from random import choice, sample
from lib.curves import create_curve_with_light_bend_and_noise
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
        end_cell = sub_grid.get_cell_at_index(0, 5)
        start_cell = sub_grid.get_cell_at_index(5, 9)

        vsk.point(start_cell.center.x, start_cell.center.y)
        vsk.point(end_cell.center.x, end_cell.center.y)

        draw_grid(vsk, sub_grid)

        unit_trajectory = Vector.from_two_points(start_cell.center, end_cell.center).normalize()
        length = Vector.from_two_points(start_cell.center, end_cell.center).length

        vsk.stroke(2)
        self.create_branched_curves(vsk, start_cell.center, length, unit_trajectory, 0, 0.7) 

    def create_branched_curves(
        self,
        vsk: vsketch.Vsketch,
        start_point: Point,
        length: float,
        unit_trajectory: Vector,
        count: int,
        trunk_width: float,
    ):
        if count > 2:
            return

        # Draw the curve as described by the given parameters 
        bend_clockwise = unit_trajectory.x > 0
        end_point = add(start_point, unit_trajectory * length)
        curve = create_curve_with_light_bend_and_noise(
            vsk,
            start_end = (start_point, end_point),
            control_point_locations = (0.25, 0.75),
            trunk_width = trunk_width,
            tip_width = 0.0,
            control_point_precision = 0.2,
            bend_clockwise=bend_clockwise,
        )
        vsk.polygon(curve[0], curve[1])

        # Pick N random points on the curve to start branching from
        # and recurse
        indices_first_fourth = sample(range(0, len(curve[0]) // 4), 1)
        indices_second_fourth = sample(range(len(curve[0]) // 4, len(curve[0]) // 2), 1)
        indices_third_fourth = sample(range(len(curve[0]) // 2, 3 * len(curve[0]) // 4), 1)
        indices_fourth_fourth = sample(range(3 * len(curve[0]) // 4, len(curve[0]) - 1), 1)
        indices = indices_first_fourth + indices_second_fourth + indices_third_fourth + indices_fourth_fourth
        start_points = [curve[0][indices], curve[1][indices]]
        # draw_vector(vsk, unit_trajectory * length, start_point)
        for i in range(len(start_points[0])):
            # Pick a random angle to rotate
            angle = choice([30, 20])
            new_length = length * 0.6
            new_count = count + 1
            rotated_trajectory = rotate_vector(unit_trajectory, angle, degrees = True)
            self.create_branched_curves(
                vsk,
                Point(start_points[0][i], start_points[1][i]),
                new_length,
                rotated_trajectory,
                new_count,
                trunk_width * 0.50
            )
        
        # TODO: Instead of just branching out from one point, select N random points 
        # along the curve and branch out from there.

        # TODO: Fix the width of the branches such that the branching feels a bit more natural.

