from math import copysign
from random import uniform, choice
import sys
import bezier
from lib.curves import create_curve_with_light_bend, create_curve_with_light_bend_and_noise, draw_cubic_bezier
from lib.point_utils import add
from lib.random_utils import random_n_elements_across_k_partitions
import vsketch
import numpy as np
from typing import Optional, List 
from lib.grid import Grid, create_grid_with_padding, draw_grid
from lib.vector import Vector, draw_vector, rotate_vector
from shapely.geometry import LineString, Point

# Next Steps
# - Use the builder pattern for constructing the leaf (width, length, orientation, num_veins, vein strategy, sub-vein strategy, etc.)
# - Recursively generate sub-veins: https://natureofcode.com/book/chapter-8-fractals/
# - Add randomness to the shape of the center vein, the direction of the veins, the outer shape of the leaf
# - I might want to use the intersection function to actually find the point on the curve. I might just be able to use `evaluate_multi` though

class LeafSketch(vsketch.SketchClass):
    leaf_length = vsketch.Param(3.0)
    leaf_width = vsketch.Param(2.0)
    num_veins = vsketch.Param(5)

    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("9in", "11in", landscape=True, center=False)
        vsk.scale("in")

        grid = create_grid_with_padding(0.5, (0.50, 0.50)) 
        # draw_grid(vsk, grid)
        
        leaf_base = grid.get_cell_at_index(5, 2).center
        leaf = Leaf(leaf_base, self.leaf_length, self.leaf_width, self.num_veins, vsk)
        draw_leaf(vsk, leaf)

        leaf_base = grid.get_cell_at_index(5, 8).center
        leaf = Leaf(leaf_base, self.leaf_length, self.leaf_width, self.num_veins, vsk)
        draw_leaf(vsk, leaf)
    
        leaf_base = grid.get_cell_at_index(14, 8).center
        leaf = Leaf(leaf_base, self.leaf_length, self.leaf_width, self.num_veins, vsk)
        draw_leaf(vsk, leaf, 3)

        leaf_base = grid.get_cell_at_index(14, 14).center
        leaf = Leaf(leaf_base, self.leaf_length, self.leaf_width, self.num_veins, vsk)
        draw_leaf(vsk, leaf, 3)

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")



class Leaf():
    """
    Leaf is an abstraction for drawing leaves. 

    A `vsketch.Vsketch` object can be passed in for debugging

    Limitations:
    - Currently the leaf can only be drawn in a portrait / vertical orientation. This class needs to bezier
      updated in order to support arbitrary orientations. The internal Bezier curves will have to be constructed
      with orientation in mind. Perhaps there is an additional rotation step that can take place to allow for
      different orientations (instead of doing it at the layer of Bezier construction)
    - The y position of the control points are currently hard coded
    """
    
    # The veins will be bounded by another bezier curve that is structurally equivalent to the
    # main leaf curves, except that the control point influenced by width will be adjusted
    # down by this fraction
    VEIN_BOUNDARY_PERCENTAGE_OF_WIDTH = 0.9

    # The trajectory of the veins
    # TODO: Add some randomness to this
    VEIN_TRAJECTORY = Vector(-1.0, -1.0).normalize()

    def __init__(self, base: Point, length: float, width: float, num_veins: int = 5, vsk: Optional[vsketch.Vsketch] = None):
        if vsk:
            self.vsk = vsk

        self.base = base
        self.tip = Point(base.x, base.y - length)   
        self.width = width
        self.length = length
        self.num_veins = num_veins

        # right_vein_trajectory = Vector(1.0, -1.0).normalize() 
        self.right_side = self._construct_side(self.width, left = False)
        self.left_side = self._construct_side(self.width)

        trajectory = Vector.from_two_points(self.base, self.tip).normalize()

        self.veins = self.create_branched_curves(
            self.vsk,
            self.base,
            length,
            trajectory,
            count = 0,
            trunk_width = 0.2,
            control_point_precision = 0.01,
        )

        # self.right_veins = self._construct_veins(right_vein_trajectory)
        # left_vein_trajectory = Vector(-1.0, -1.0).normalize() 
        # self.left_veins = self._construct_veins(left_vein_trajectory)


    def _construct_side(self, width: float, left: bool = True) -> bezier.Curve:
        """
        Constructs a side of the leaf. direction == True will draw the left side. False draws the right
        """
        modifier = -1 if left else 1 
        width *= modifier
        nodes = np.asarray([
            [self.base.x, self.base.x + width, self.base.x, self.tip.x],
            [self.base.y, self.base.y - 1.0, self.tip.y + 1.0, self.tip.y]
        ])
        return bezier.Curve(nodes, degree = 3)


    def create_branched_curves(
        self,
        vsk: vsketch.Vsketch,
        start_point: Point,
        length: float,
        unit_trajectory: Vector,
        count: int,
        trunk_width: float,
        control_point_precision: float,
    ) -> List[bezier.Curve]:
        curves = []
        if count > 2:
            return [] 

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
        curves.append(no_noise_curve)
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
            recursed_curves = self.create_branched_curves(
                vsk,
                Point(start_points[0][i], start_points[1][i]),
                new_length,
                rotated_trajectory,
                new_count,
                trunk_width * 0.6,
                control_point_precision * 0.6,
            )
            curves = curves + recursed_curves
        return curves

   
def draw_leaf(vsk: vsketch.Vsketch, leaf: Leaf, weight: int = 1):
    vsk.strokeWeight(weight)

    draw_cubic_bezier(vsk, leaf.left_side)
    draw_cubic_bezier(vsk, leaf.right_side)

    vsk.stroke(1)

    for vein in leaf.veins:
        draw_cubic_bezier(vsk, vein)
        
    # combined = leaf.left_veins + leaf.right_veins
    # for vein in combined:
    #     if vein.level == 1:
    #         vsk.strokeWeight(max(1, weight - 2))
    #     else: 
    #         vsk.strokeWeight(weight)
    #     draw_cubic_bezier(vsk, vein.curve)

    vsk.strokeWeight(weight)

    # vsk.line(
    #     leaf.base.x, leaf.base.y, leaf.base.x, leaf.base.y - leaf.length)


if __name__ == "__main__":
    LeafSketch.display()
