import bezier
import vsketch
import numpy as np
from typing import Sequence, Tuple
from lib.grid import create_grid_with_padding
from lib.vector import Vector, VectorTextOptions, VectorTextOrientation, draw_vector
from shapely.geometry import Point


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
        vsk.size("letter", landscape=True, center=False)
        vsk.scale("in")

        grid = create_grid_with_padding(0.5, (0.25, 0.25)) 
        # draw_grid(vsk, grid)
        center = grid.center
    
        leaf_base = center
        leaf_tip = Point(center.x, center.y - self.leaf_length)   

        leaf = Leaf(leaf_base, self.leaf_length, self.leaf_width, self.num_veins)
        draw_leaf(vsk, leaf)

        # TODO: Introduce some randomness here
        vsk.line(
            leaf_base.x, leaf_base.y, leaf_tip.x, leaf_tip.y)


    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")


class Leaf():
    """
    Leaf is an abstraction for drawing leaves. 

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

    def __init__(self, base: Point, length: float, width: float, num_veins: int = 5):
        self.base = base
        self.tip = Point(base.x, base.y - length)   
        self.width = width
        self.length = length
        self.num_veins = num_veins
        self.left_side = self._construct_side(-1 * self.width)
        self.right_side = self._construct_side(self.width)
        self.left_vein_boundary = self._construct_side(-1 * self.VEIN_BOUNDARY_PERCENTAGE_OF_WIDTH * self.width)
        self.left_veins = self._construct_veins()


    def _construct_side(self, width: float) -> bezier.Curve:
        nodes = np.asarray([
            [self.base.x, self.base.x + width, self.base.x, self.tip.x],
            [self.base.y, self.base.y - 1.0, self.tip.y + 1.0, self.tip.y]
        ])
        return bezier.Curve(nodes, degree = 3)

    
    def _construct_veins(self) -> Sequence[Tuple[Point, Point]]:
        vein_gap = (self.base.y - self.tip.y) / self.num_veins
        veins = []
        for i in range(self.num_veins):
            vein_start = Point(self.base.x, self.base.y - vein_gap * i)

            # Project vein vector onto the leaf width to get an appropriately sized
            # vein (i.e fits in the leaf)
            vein_trajectory = Vector(-1.0, -1.0).get_unit_vector() 
            vein_end = self._compute_vein_end(vein_start, vein_trajectory, self.left_vein_boundary)
            veins.append((vein_start, vein_end))
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
    

def draw_leaf(vsk: vsketch.Vsketch, leaf: Leaf):
    left_nodes = leaf.left_side.nodes
    right_nodes = leaf.right_side.nodes
    assert(left_nodes.shape == (2, 4)) # We expect a cubic Bezier curve
    assert(right_nodes.shape == (2, 4)) # We expect a cubic Bezier curve

    def sketch_bezier(nodes: np.ndarray):
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
    
    sketch_bezier(left_nodes)
    sketch_bezier(right_nodes)
    
    vsk.stroke(1)
    # DEBUG
    # sketch_bezier(leaf.left_vein_boundary.nodes)
        
    for vein in leaf.left_veins:
        vein_start = vein[0]
        vein_end = vein[1]
        vsk.line(vein_start.x, vein_start.y, vein_end.x, vein_end.y)



if __name__ == "__main__":
    LeafSketch.display()
