from __future__ import annotations

import vsketch
import numpy as np
from lib.vector import Vector, draw_vector
from lib.grid import create_grid_with_padding, draw_grid
from lib.point_utils import add
from typing import List
from shapely import Point


# 11/13
# - Implemented a branching algorithm that selects arbitrary nodes in the existing set and starts branching at an angle perpendicular
#   to where we are at now
# - We're not getting the space filling that we want. We're getting a lot of empty space. Hopefully this will be fixed when the collition
#   detection algorithm is in place.
# # TODO: Vary the perturbation more as the branch level increases
# # TODO: Vary the decrease in the radius more as the branch level increases
# # TODO: Implement collision detection
# # TODO: There may be a bug where when you create a perpendicular branch, the start of the branch is not from the center of the previous node
#         by design and so you see a small space when we branch. This might be fixed if we select the right starting radius of the circles to be small enough
#

# 10/29
# - Currently the space filling algorithm initializes nodes and then creates new nodes from those initial nodes
#   by taking the direction of travel, perturbating it, and creating new nodes on the perimeter untill we reach the boundary
# Next steps
# - Create a collision detection algorithm to ensure no new nodes overlap with any existing ones (attempt to implement this with matrix multiplication / vectorization)
# - Implement a branching algorithm (reducing the width, increasing the pertubation range)
# - Ensure that we are space filling well


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

    draw_circles = vsketch.Param(True)
    draw_vector_tips = vsketch.Param(True)

    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("9in", "11in", landscape=True, center=False)
        vsk.scale("in")
        grid = create_grid_with_padding(0.5, (0.50, 0.50))
        center = grid.center
        max_circle_radius = min(abs(grid.right - center.x), abs(grid.top - center.y))

        bfill = BranchFill(Circle(center.x, center.y, max_circle_radius))
        bfill.fill()
        self.draw_bfill(vsk, bfill)

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")

    def draw_bfill(self, vsk: vsketch.Vsketch, bfill: BranchFill):
        boundary = bfill.boundary
        vsk.circle(boundary.x, boundary.y, boundary.r * 2.0)

        # for node in bfill.initial_nodes:
        #     if self.draw_circles:
        #         vsk.circle(node.circle.x, node.circle.y, node.circle.r * 2.0)
        #     self.draw_node(vsk, node)

        for node in bfill.all_nodes:
            if self.draw_circles:
                vsk.circle(node.circle.x, node.circle.y, node.circle.r * 2.0)
            self.draw_node(vsk, node)

    def draw_node(self, vsk: vsketch.Vsketch, node: Node):
        if self.draw_circles:
            vsk.circle(node.circle.x, node.circle.y, node.circle.r * 2.0)

        # Draw a single line segment in the direction of node.direction, passing
        # through the center of the circle, of length equal to the diameter
        # of the circle
        reverse_direction = node.direction.normalize() * -1.0 * node.circle.r
        start_point = add(node.circle.center, reverse_direction)
        draw_vector(vsk, node.direction * 2.0, start_point, draw_tip=self.draw_vector_tips)


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
    def __init__(self, circle: Circle, direction: Vector, branch_level: int = 0):
        self.circle = circle
        self.direction = direction
        self.branch_level = branch_level

    def is_within(self, other: Circle) -> bool:
        return self.circle.is_within(other)

    @staticmethod
    def from_previous_node(
            previous_node: Node,
            new_direction: Vector,
            scale_radius_down_by: float = 1.0,
            branch_level: int = 0,
    ) -> Node:
        """
        Creates a new node from the previous node by pertubating the direction of the previous node
        and scaling down the radius of the new node.
        @param previous_node: The previous node to use as a reference to create the new node. We primarily use the previous
        node's circle's center and radius to create the new node.
        @param new_direction: The direction in which to create the new node.
        @param scale_radius_down_by: How much to scale down the radius of the new node by. This is a value between 0 and 1.
        @param start_from_perimeter: If True, then the new node's will begin from the perimeter of the previous node. If False,
        then the new node will begin from previous node's center.
        @param is_branch: If True, then the new node will have a branch_level of previous_node.branch_level + 1. If False,
        """

        # Using the current nodes center and direction, perturbate the angle a little bit
        # and find a point on the perimeter of the circle. This will be the tangent point
        # where the new node's circle will intersect with the current node.
        #
        # First perturbate to find the tangent point. We normalize the direction and scale it by the previous
        # nodes radius to ensure we're on the perimeter of the circle
        tangent_point = add(previous_node.circle.center, new_direction.normalize() * previous_node.circle.r)
        starting_point = tangent_point

        # Second, scale down the new radius, and find the new center
        radius = previous_node.circle.r * scale_radius_down_by
        scaled_direction = new_direction.normalize() * radius
        new_center = add(starting_point, scaled_direction)
        circle = Circle(new_center.x, new_center.y, radius)

        return Node(circle, scaled_direction, branch_level)

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
        self.all_nodes: List[Node] = []
        self.NUM_INITIAL_NODES = 3
        self.INITIAL_NODE_RADIUS_PCT = 0.01
        self.NEXT_NODE_RADIUS_REDUCTION_PCT = 0.99
        self.NEXT_NODE_CROSS_BOUNDARY_RETRIES = 3
        self.MIN_NODE_RADIUS_THRESHOLD_INCHES = 0.01

    def fill(self):
        self.initial_nodes = self.select_initial_nodes()
        for node in self.initial_nodes:
            initial_pertubation = (np.deg2rad(-15.0), np.deg2rad(15.0))
            self.all_nodes.extend(self.begin_fill_from_node(node, node.direction, initial_pertubation))
            self.all_nodes.append(node)

        # Branching algorithm:
        # Randomly pick nodes from self.all_nodes, pick an angle that's perpendicular (with some perturbation)
        # to the direction of the original node's travel, and then continue filling but with a tighter radius,
        # and a larger perturbation range.
        NUM_BRANCH_STARTS = 80
        for i in range(NUM_BRANCH_STARTS):
            # Randomly pick a node, and a direction that's perpendicular to the node's direction. Perturbate the direction
            # a bit
            node = np.random.choice(self.all_nodes)
            # Chose the direction to branch and create the branching node
            clockwise = np.random.choice([True, False])
            direction_to_branch = node.direction.normalize().get_perpendicular(clockwise=clockwise)

            # Begin a fill from this node
            self.all_nodes.extend(self.begin_fill_from_node(node, direction_to_branch, initial_pertubation))




    def begin_fill_from_node(self, node: Node, initial_direction: Vector, pertubation_range: (float, float)) -> List[Node]:
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
        next_branch_level = node.branch_level + 1
        current_node = node
        current_direction = initial_direction
        retries = 0
        while True:
            # Using the current nodes center and direction, perturbate the angle a little bit
            # and find a point on the perimeter of the circle. This will be the tangent point
            # where the new node's circle will intersect with the current node.
            #
            # First perturbate to find the tangent point
            direction_pertubated = current_direction.perturbate(pertubation_range)
            new_node = Node.from_previous_node(current_node, direction_pertubated, scale_radius_down_by=self.NEXT_NODE_RADIUS_REDUCTION_PCT, branch_level=next_branch_level)

            # If the new circle is out of the boundary, allow for some retries with pertubations, but
            # eventually give up and return.
            if not new_node.is_within(self.boundary):
                if retries > self.NEXT_NODE_CROSS_BOUNDARY_RETRIES:
                    break
                retries += 1
                continue

            # If the new node is too small, then give up
            if new_node.circle.r < self.MIN_NODE_RADIUS_THRESHOLD_INCHES:
                break

            # Register the new node
            new_nodes.append(new_node)
            current_node = new_node
            current_direction = direction_pertubated
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
