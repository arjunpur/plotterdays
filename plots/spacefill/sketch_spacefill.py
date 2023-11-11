from __future__ import annotations

import vsketch
import numpy as np
from lib.vector import Vector, draw_vector
from lib.grid import create_grid_with_padding, draw_grid
from lib.point_utils import add
from typing import List
from shapely import Point


# 10/29
# - Currently the space filling algorithm initializes nodes and then creates new nodes from those initial nodes
#   by taking the direction of travel, perturbating it, and creating new nodes on the perimeter untill we reach the boundary
# Next steps
# - Create a collision detection algorithm to ensure no new nodes overlap with any existing ones (attempt to implement this with matrix multiplication / vectorization)
# - Implement a branching algorithm (reducing the width, increasing the pertubation range)
# - Ensure that we are space filling well
# TODO: There's a bug where the next_node's aren't completely non-intersecting (THE MATH IS WRONG!!!!!!!!!!!!!!)


# Space Branch-Fill Algorithm
# - Start with some initial nodes within a circle. These nodes are circles with a *radius*, and a *direction of travel*.
# - We are then going to pick a point on the perimeter of the circle using the *direction of travel* and a
#   slight angle pertubation (proportional to the *branch level / width* - thicker branches have less sway to them)
# - Continue doing this untill you get to a point where you cannot create a new circle without overlapping with an existing circle / boundary (with 3 retries of pertubations).
# - Now to branch, pick a random node and repeat the above process (with a smaller branch width - thus more pertubations)
# - Keep doing this until most random selections of nodes fail to create a new circle (with 3 retries of pertubations)


# Collision Detection
# Q: How do we detect collisions between circles?
# Q: How do we even represent the circles?
#  - I might just be able to represent them as (x, y, radius) tuples instead of using a library like shapely
#  - I can then use the distance formula to check if the distance between the centers of the circles is less than the sum of the radii

class SpaceFillSketch(vsketch.SketchClass):
    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("9in", "11in", landscape=True, center=False)
        vsk.scale("in")
        grid = create_grid_with_padding(0.5, (0.50, 0.50))
        center = grid.center
        max_circle_radius = min(abs(grid.right - center.x), abs(grid.top - center.y))

        bfill = BranchFill(Circle(center.x, center.y, max_circle_radius))
        bfill.fill()
        draw(vsk, bfill)

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")


def draw(vsk: vsketch.Vsketch, bfill: BranchFill):
    boundary = bfill.boundary
    vsk.circle(boundary.x, boundary.y, boundary.r * 2.0)

    for node in bfill.initial_nodes:
        vsk.circle(node.circle.x, node.circle.y, node.circle.r * 2.0)
        draw_vector(vsk, node.direction, node.circle.center, draw_tip=True)

    for node in bfill.all_nodes:
        vsk.circle(node.circle.x, node.circle.y, node.circle.r * 2.0)
        draw_vector(vsk, node.direction, node.circle.center, draw_tip=True)


class Circle:
    def __init__(self, x, y, r):
        self.x = x
        self.y = y
        self.r = r
        self.np = np.array([x, y])

    def is_within(self, other: Circle) -> bool:
        d = np.linalg.norm(self.np - other.np)
        # Check if circle2 is within circle1
        return d + self.r <= other.r

    def intersects(self, other: Circle) -> bool:
        d = np.linalg.norm(self.np - other.np)
        # Compare distance with sum of radii
        return d < (self.r + other.r)

    @property
    def center(self) -> Point:
        return Point(self.x, self.y)

    def random_circles_within(self, radius: float, n: int = 3) -> List[Circle]:
        """Generate n random, non-overlapping circles within this circle"""

        attempts = 0
        # TODO: Ensure that these points are non-overlapping.
        circles = []
        while len(circles) < n:
            attempts += 1
            if attempts >= 20:
                raise Exception("Could not generate non-overlapping circles within the given number of attempts")

            # Generate a random angle and distance from the center
            angle = 2 * np.pi * np.random.rand()
            # To ensure a uniform distribution of points inside the circle, we use np.sqrt(np.random.rand()).
            # The reason is that the area of a circle grows with the square of its radius.
            # By taking the square root of a uniformly distributed random number, we are compensating for
            # that squared growth, resulting in an even distribution of points.
            distance = self.r * np.sqrt(np.random.rand())

            # Convert polar coordinates to cartesian coordinates
            x = self.x + distance * np.cos(angle)
            y = self.y + distance * np.sin(angle)
            circle = Circle(x, y, radius)

            intersects_with_existing = False
            for existing_circle in circles:
                if circle.intersects(existing_circle):
                    intersects_with_existing = True

            if not intersects_with_existing:
                circles.append(circle)

        return circles


class Node:
    def __init__(self, circle: Circle, direction: Vector):
        self.circle = circle
        self.direction = direction


class AllCircles:
    def __init__(self):
        # Represent circles as numpy array of (x, y, r) tuples and then write the intersection function in
        # terms of numpy array operations
        self.circles = np.array([])

    def does_circle_intersect(self, circle: Circle) -> bool:
        pass


class BranchFill:
    def __init__(self, boundary: Circle):
        # TODO: Make this an abstract class that has the key methods needed for the boundary (ex. is_within)
        # This will allow us to use different shapes as the boundary.
        self.boundary = boundary
        self.initial_nodes = None
        self.all_nodes = []
        self.NUM_INITIAL_NODES = 3
        self.INITIAL_NODE_RADIUS_PCT = 0.1
        self.NEXT_NODE_RADIUS_REDUCTION_PCT = 0.95
        self.NEXT_NODE_CROSS_BOUNDARY_RETRIES = 3

    def fill(self):
        self.initial_nodes = self.select_initial_nodes()
        for node in self.initial_nodes:
            initial_pertubation = (np.deg2rad(-15.0), np.deg2rad(15.0))
            self.all_nodes.extend(self.begin_fill_from_node(node, initial_pertubation))

    def begin_fill_from_node(self, node: Node, pertubation_range: (float, float)) -> List[Node]:
        """
        Starts from a node's center and uses the direction to repeatedly create new nodes till it
        reaches the boundary and cannot create more nodes (with 3 retries of pertubations).
        Each new node is created by pertubating the angle of the direction vector and finding a point
        on the perimeter of the node's circle. This point is then used to create a new circle with a
        smaller radius such that the point is a tangent to the older node.
        @param node:
        @param pertubation_range:
        @return:
        """
        new_nodes = []
        current_node = node
        retries = 0
        while True:
            # Using the current nodes center and direction, perturbate the angle a little bit
            # and find a point on the perimeter of the circle. This will be the tangent point
            # where the new node's circle will intersect with the current node.
            #
            # First perturbate to find the tangent point
            direction_perturbated = current_node.direction.perturbate(pertubation_range)
            tangent_point = add(current_node.circle.center, direction_perturbated)

            # Second, scale down the new radius, and find the new center
            radius = current_node.circle.r * self.NEXT_NODE_RADIUS_REDUCTION_PCT
            new_direction = direction_perturbated.normalize() * radius
            center = add(tangent_point, new_direction)
            circle = Circle(center.x, center.y, radius)

            # If the new circle is out of the boundary, allow for some retries with pertubations, but
            # eventually give up and return.
            if not circle.is_within(self.boundary):
                if retries > self.NEXT_NODE_CROSS_BOUNDARY_RETRIES:
                    break
                retries += 1
                continue

            # Register the new node
            new_node = Node(circle, new_direction)
            new_nodes.append(new_node)
            current_node = new_node
        return new_nodes

    def select_initial_nodes(self) -> List[Node]:
        """
        Selects the initial nodes from an inner circle defined by `self.INITIAL_NODE_RADIUS_PCT * self.boundary.r`. Randomly selects
        the points in the circles, and directions for the nodes to travel in
        """
        inner_circle = Circle(self.boundary.x, self.boundary.y, self.boundary.r * 0.4)
        initial_node_radius = self.boundary.r * self.INITIAL_NODE_RADIUS_PCT

        random_circles_within = inner_circle.random_circles_within(initial_node_radius, self.NUM_INITIAL_NODES)

        random_angles = [2 * np.pi * np.random.rand() for _ in range(len(random_circles_within))]
        initial_nodes = []
        for (i, point) in enumerate(random_circles_within):
            circle = Circle(point.x, point.y, initial_node_radius)
            direction = Vector.from_length_and_angle(initial_node_radius, random_angles[i])
            node = Node(circle, direction)
            initial_nodes.append(node)
        return initial_nodes

    def within_boundary(self, circle: Circle) -> bool:
        pass


# The goal is to create an algorithm that fills a space with non-overlapping circles by creating circles
# on the perimeter of each created circle in a direction of travel. As you travel further, the radius of the
# newly created circle decreases. The angle in which you create the circle is also pertubed. As the radius
# gets smaller, the angle pertubations also increase.

# To get branches, you pick a node at random and try to create a new node


# Define the boundary of the space to fill with branches
# - This will likely be a circle. Also define an inner circle where the initial nodes will be placed.

# Randomly select the initial nodes from the inner circle.
# For each initial node, start the branch-fill algorithm given some general direction and other initial conditions like
# - Starting radius

# For each node, create line segments with random angle pertubations till the line segment intersects with the boundary of the radius defined.

# Now create branches by:
#


if __name__ == "__main__":
    SpaceFillSketch.display()
