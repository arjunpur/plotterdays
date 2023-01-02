from math import copysign
from random import uniform, choice
import sys
import bezier
from lib.curves import create_curve_with_light_bend, create_curve_with_light_bend_and_noise, draw_cubic_bezier
from lib.point_utils import add
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

        
        # num_lines = 250
        # x_coords = np.linspace(0., 250., 1000)
        # perlin = vsk.noise(x_coords * 0.1, np.linspace(0, 5., num_lines))
        # print(perlin.shape)
        # for j in range(num_lines):
        #     vsk.polygon(x_coords, j + perlin[:, j] * 10)


    def draw_branch(self, grid: Grid, vsk: vsketch.Vsketch) -> None:
        branch_start = grid.get_cell_at_index(10, 3).center
        branch_end = grid.get_cell_at_index(10, 18).center
        
        num_leaves = 5
        leaf_distance = (branch_end.x - branch_start.x) / num_leaves

        vsk.line(branch_start.x, branch_start.y, branch_end.x, branch_end.y)
        
        for i in range(1, num_leaves):
            leaf_base = Point(branch_start.x + leaf_distance * i, branch_start.y)
            leaf = Leaf(leaf_base, self.leaf_length, self.leaf_width, self.num_veins, vsk)
            draw_leaf(vsk, leaf)

        vsk.line(branch_start.x, branch_start.y, branch_end.x, branch_end.y)


    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")


class Vein():

    def __init__(
        self,
        vein_start: Point,
        vein_trajectory: Vector,
        vein_boundaries: List[bezier.Curve],
        level: int = 0,
        vsk: Optional[vsketch.Vsketch] = None,
    ):
        self.level = level
        if vsk:
            self.vsk = vsk

        for boundary in vein_boundaries:
            computed_vein_end = self._compute_vein_end(vein_start, vein_trajectory, boundary)
            if computed_vein_end:
                self.vein_end = computed_vein_end
                break
        if not self.vein_end:
            raise Exception("could not compute vein end - invalid boundaries")

        vein_line = LineString([vein_start, self.vein_end])
        self.vein_boundary = vein_boundaries
        
        # TODO: This might not work if the orientation of the vein is reverse
        trajectory = Vector.from_two_points(vein_start, self.vein_end).normalize()
        left_or_right = trajectory.x < 0
        self.curve = create_curve_with_light_bend(
            (vein_start, self.vein_end),
            (0.25, 0.75),
            vein_line.length * 0.2,
            vein_line.length * 0.1,
            vein_line.length * 0.1,
            bend_clockwise = left_or_right,
        )
        

    def _compute_vein_end(self, vein_start: Point, vein_trajectory: Vector, side_curve: bezier.Curve) -> Optional[Point]:
        try:
            scaled_vein_trajectory = vein_trajectory * sys.maxsize
            bounding_point_end = add(vein_start, scaled_vein_trajectory)
            nodes = np.asarray([[vein_start.x, bounding_point_end.x], [vein_start.y, bounding_point_end.y]])
            line_segment = bezier.Curve(
                nodes,
                degree = 1
            )
            intersection_s_value_with_vein_boundary = side_curve.intersect(line_segment)[0][0]
            intersection_point = side_curve.evaluate(intersection_s_value_with_vein_boundary)
            return Point(intersection_point)
        except:
            return None

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

        right_vein_trajectory = Vector(1.0, -1.0).normalize() 
        self.right_side = self._construct_side(self.width, left = False)
        self.right_veins = self._construct_veins(right_vein_trajectory)

        left_vein_trajectory = Vector(-1.0, -1.0).normalize() 
        self.left_side = self._construct_side(self.width)
        self.left_veins = self._construct_veins(left_vein_trajectory)


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

   
    def _construct_veins(self, vein_trajectory: Vector) -> List[Vein]:
        """
        Constructs a series of veins for a particular side of the leaf. Each vein is modelled as a 
        Bezier curve.
        """
        vein_gap = abs(self.base.y - self.tip.y) / self.num_veins
        veins = []
        
        left_or_right = vein_trajectory.x < 0
        vein_boundary = self._construct_side(self.VEIN_BOUNDARY_PERCENTAGE_OF_WIDTH * self.width, left = left_or_right)
        last_vein = None
        for i in range(1, self.num_veins):
            vein_start = Point(self.base.x, self.base.y - vein_gap * i)
            vein = Vein(vein_start, vein_trajectory, [vein_boundary])
            veins.append(vein)

            # Construct sub-veins at hard coded points on the main vein. Rotate the main vein'scaled_vein_trajectory
            # trajectory by angles defined below
            # The sub-vein intersects either with the main vein boundary, or the last main vein created. 
            # TODO: Replace with better collision detection
            parameters = np.asarray([
                uniform(0.1, 0.3),
                uniform(0.35, 0.5),
                uniform(0.55, 0.70),
                uniform(0.75, 0.85),
                uniform(0.75, 0.95),
            ])
            points = vein.curve.evaluate_multi(parameters)
            rotations_for_trajectories = np.asarray([
                uniform(-30.0, 0.0),
                uniform(-20.0, 30.0),
                uniform(20.0, 50.0),
                uniform(20.0, 50.0),
            ]) * copysign(1, vein_trajectory.x) 
    
            rotator = np.vectorize(lambda angle: rotate_vector(vein_trajectory, angle, degrees = True))

            for j, trajectory in enumerate(rotator(rotations_for_trajectories)):
                sub_vein_start = Point(points[0][j], points[1][j])
                boundaries = [vein_boundary]
                if last_vein:
                    boundaries = [last_vein.curve, vein_boundary]
                sub_vein = Vein(sub_vein_start, trajectory, boundaries, level = 1, vsk = self.vsk)
                veins.append(sub_vein)

            last_vein = vein

            # TODO: Iterate through these sub_vein starting points and re-call the light bend curve

            # TODO: Rotate the vein_trajectory to aim the sub_vein + add some randomness

        return veins

def draw_leaf(vsk: vsketch.Vsketch, leaf: Leaf, weight: int = 1):
    vsk.strokeWeight(weight)

    draw_cubic_bezier(vsk, leaf.left_side)
    draw_cubic_bezier(vsk, leaf.right_side)

    vsk.stroke(1)
        
    combined = leaf.left_veins + leaf.right_veins
    for vein in combined:
        if vein.level == 1:
            vsk.strokeWeight(max(1, weight - 2))
        else: 
            vsk.strokeWeight(weight)
        draw_cubic_bezier(vsk, vein.curve)

    vsk.strokeWeight(weight)

    vsk.line(
        leaf.base.x, leaf.base.y, leaf.base.x, leaf.base.y - leaf.length)



if __name__ == "__main__":
    LeafSketch.display()
