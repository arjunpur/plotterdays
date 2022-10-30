import random
import bezier
import vsketch
import numpy as np
from typing import Sequence, Optional 
from lib.grid import create_grid_with_padding, draw_grid
from lib.vector import Vector
from shapely.geometry import LineString, Point, Polygon


# Next Steps
# - Abstract away a Bezier curve and a leaf using this package: https://bezier.readthedocs.io/en/stable/python/reference/bezier.curve.html
# - Using the abstraction of the Bezier curve, use the `evaluate_multi(s)` function to get the corresponding
#   points on the curve. Use this to determine the distance between the curve and the center vein of the leaf.
#   That value will be used to scale the leaf vein.
# - Recursively generate sub-veins
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
        
        # leaf_base = grid.get_cell_at_index(5, 2).center
        # leaf = Leaf(leaf_base, self.leaf_length, self.leaf_width, self.num_veins, vsk)
        # draw_leaf(vsk, leaf)


        # leaf_base = grid.get_cell_at_index(5, 8).center
        # leaf = Leaf(leaf_base, self.leaf_length, self.leaf_width, self.num_veins, vsk)
        # draw_leaf(vsk, leaf)


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


        # TODO: Introduce some randomness here


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
    VEIN_BOUNDARY_PERCENTAGE_OF_WIDTH = 0.7

    # The trajectory of the veins
    # TODO: Add some randomness to this
    VEIN_TRAJECTORY = Vector(-1.0, -1.0).get_unit_vector()

    def __init__(self, base: Point, length: float, width: float, num_veins: int = 5, vsk: Optional[vsketch.Vsketch] = None):
        if vsk:
            self.vsk = vsk

        self.base = base
        self.tip = Point(base.x, base.y - length)   
        self.width = width
        self.length = length
        self.num_veins = num_veins
        self.right_side = self._construct_side(self.width, left = False)
        self.right_vein_boundary = self._construct_side(self.VEIN_BOUNDARY_PERCENTAGE_OF_WIDTH * self.width, left = False)
        self.right_veins = self._construct_veins(self.right_vein_boundary, left = False)

        self.left_side = self._construct_side(self.width)
        self.left_vein_boundary = self._construct_side(self.VEIN_BOUNDARY_PERCENTAGE_OF_WIDTH * self.width)
        self.left_veins = self._construct_veins(self.left_vein_boundary)


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

   
    def _construct_veins(self, vein_boundary: bezier.Curve, left: bool = True) -> Sequence[bezier.Curve]:
        """
        Constructs a series of veins for a particular side of the leaf. Each vein is modelled as a 
        Bezier curve.
        """
        vein_gap = (self.base.y - self.tip.y) / self.num_veins
        veins = []

        # Project vein vector onto the leaf width to get an appropriately sized
        # vein (i.e fits in the leaf)
        modifier = -1 if left else 1 
        vein_trajectory = Vector(1.0 * modifier, -1.0).get_unit_vector() 

        for i in range(1, self.num_veins):
            vein_start = Point(self.base.x, self.base.y - vein_gap * i)
            vein_end = self._compute_vein_end(vein_start, vein_trajectory, vein_boundary)

            # Pick two random points within a semi-circle using the start / end as the 
            # ends of the diameter of the semi-circle
            # https://shapely.readthedocs.io/en/stable/manual.html#object.buffer
            # A negative number for buffer will choose the left side, whereas positive chooses the right
            vein_line = LineString([vein_start, vein_end])
            buffer_length = vein_line.length * 0.4
            vein_buffer = Polygon(LineString([vein_start, vein_end]).buffer(buffer_length * modifier, single_sided = True))

            control_point_1 = generate_random_points(1, vein_buffer)[0] 
            control_point_2 = generate_random_points(1, vein_buffer)[0] 

            nodes = np.asarray([
                [vein_start.x, control_point_1.x, control_point_2.x, vein_end.x],
                [vein_start.y, control_point_1.y, control_point_2.y, vein_end.y]
            ])
            vein = bezier.Curve(
                nodes,
                degree = 3
            )
            veins.append(vein)
        return veins


    def _compute_vein_end(self, vein_start: Point, vein_trajectory: Vector, side_curve: bezier.Curve) -> Point:
        # Find a point definitely outside the leaf (worst case the point is the tip or base
        # and the vein_trajectory is one of the two axis
        bounding_circle_radius = max(self.width, self.length)
        scaled_vein_trajectory = vein_trajectory * bounding_circle_radius
        bounding_point_end = Point(vein_start.x + scaled_vein_trajectory.x, vein_start.y + scaled_vein_trajectory.y)
        
        nodes = np.asarray([[vein_start.x, bounding_point_end.x], [vein_start.y, bounding_point_end.y]])
        line_segment = bezier.Curve(
            nodes,
            degree = 1
        )

        intersection_parameter_with_vein_boundary = side_curve.intersect(line_segment)
        intersection_point = side_curve.evaluate(intersection_parameter_with_vein_boundary[0][0])
        return Point(intersection_point)



def generate_random_points(number: int, polygon: Polygon) -> Sequence[Point]:
    points = []
    if len(polygon.bounds) == 0:
        raise Exception("could not compute control points for vein bezier")
    minx, miny, maxx, maxy = polygon.bounds
    while len(points) < number:
        pnt = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
        if polygon.contains(pnt):
            points.append(pnt)
    return points
    

def draw_leaf(vsk: vsketch.Vsketch, leaf: Leaf):
    def sketch_cubic_bezier(curve: bezier.Curve):
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

    sketch_cubic_bezier(leaf.left_side)
    sketch_cubic_bezier(leaf.right_side)

    # sketch_cubic_bezier(leaf.left_vein_boundary)
    
    vsk.stroke(1)
    # DEBUG
    # sketch_bezier(leaf.left_vein_boundary.nodes)
        
    for vein in leaf.left_veins:
        sketch_cubic_bezier(vein)

    for vein in leaf.right_veins:
        sketch_cubic_bezier(vein)

    vsk.line(
        leaf.base.x, leaf.base.y, leaf.base.x, leaf.base.y - leaf.length)



if __name__ == "__main__":
    LeafSketch.display()
