from __future__ import annotations
import numpy as np
import vsketch
from math import cos, sin
from shapely.geometry import Point, LineString

class Cell:
    def __init__(self, position: Point, width: float, height: float):
        self.position = position
        self.width = width
        self.height = height

    # TODO: These functions can be abstracted out into a base geometry that is 
    # shared with the grid.
    @property
    def top(self) -> float:
        return self.position.y

    @property
    def bottom(self) -> float:
        return self.position.y + self.height

    @property
    def left(self) -> float:
        return self.position.x

    @property
    def right(self) -> float:
        return self.position.x + self.width

    @property
    def center(self) -> Point:
        center_x = (self.right + self.left) / 2.0
        center_y = (self.top + self.bottom) / 2.0
        return Point(center_x, center_y)


# TODO currently we proxy back and forth between numpy arrays
# I might be able to just remove this abstraction to simplify the code
# Or just use numpy arrays under the hood
class Vector:
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y

    # https://stackoverflow.com/questions/44640479/type-annotation-for-classmethod-returning-instance
    @classmethod
    def from_line_string(cls, linestring: LineString) -> Vector:
        # TODO
        return cls(0.0, 0.0)

    def to_line_string(self, start_point: Point) -> LineString:
        end_point = Point(start_point.x + self.x, start_point.y + self.y)
        return LineString([start_point, end_point])

    def to_numpy_array(self) -> np.ndarray:
        return np.array([self.x, self.y])

    @classmethod
    def from_numpy_array(cls, np_array: np.ndarray) -> Vector:
        assert(np_array.shape == (2,))
        return cls(np_array[0], np_array[1])

# TODO: Understand how vector rotation actually works mathematically
# TODO: Use shapely to perform the affine transformation. There's not
# need for me to do this with numpy directly. I'm unnecessarily using
# the Vector abstraction
def rotate_vector(
    vec: Vector,
    theta: float,
    degrees: bool = False,
) -> Vector:
    angle_radians = theta
    if degrees:
        angle_radians = np.radians(theta)

    rot_matrix = np.array([
        [cos(angle_radians), -sin(angle_radians)],
        [sin(angle_radians), cos(angle_radians)]
    ])
    rotated = np.dot(vec.to_numpy_array(), rot_matrix)
    return Vector.from_numpy_array(rotated)


def normalize(vec: Vector) -> Vector:
    ndarray = vec.to_numpy_array()
    norm = np.linalg.norm(ndarray)
    if norm == 0:
        raise Exception("could not compute norm of a 0 vector")
    return Vector.from_numpy_array(ndarray / norm)


class Grid:
    """
    A Cartesian grid that is agnostic of units. Callers are expected to maintain consistency of
    units when using the grid.

    `position` refers to the Grid's top-left coordinates.
    """

    def __init__(
        self,
        width: float,
        height: float,
        top_left: Point,
        num_cols: int,
        num_rows: int,
    ) -> None:
        self.width = width
        self.height = height
        self.top_left = top_left
        self.num_cols = num_cols
        self.num_rows = num_rows
        self._cell_width = self.width / self.num_cols
        self._cell_height = self.height / self.num_rows

    def get_cell_at_index(self, row: int, col: int) -> Cell:
        width = self._cell_width
        height = self._cell_height
        cell_x = self.top_left.x + (col * width)
        cell_y = self.top_left.y + (row * height)
        return Cell(Point(cell_x, cell_y), width, height)


    def get_vector(self, x: float, y: float) -> Vector:
        return Vector(
            x * self._cell_width,
            y * self._cell_height,
        )

    @property
    def top(self) -> float:
        return self.get_cell_at_index(0, 0).position.y

    @property
    def bottom(self) -> float:
        return self.get_cell_at_index(self.num_rows, 0).position.y

    @property
    def left(self) -> float:
        return self.get_cell_at_index(0, 0).position.x

    @property
    def right(self) -> float:
        return self.get_cell_at_index(0, self.num_cols).position.x

    @property
    def center(self) -> Point:
        center_x = (self.right + self.left) / 2.0
        center_y = (self.top + self.bottom) / 2.0
        return Point(center_x, center_y)


def draw_grid(vsk: vsketch.Vsketch, grid: Grid) -> None:
    draw_grid_border(vsk, grid)

    for row in range(1, int(grid.num_rows)):
        # Draw the column lines (vertical)
        for col in range(1, int(grid.num_cols)):
            cell = grid.get_cell_at_index(row, col)
            vsk.line(cell.position.x, grid.top, cell.position.x, grid.bottom)

        cell = grid.get_cell_at_index(row, 0)
        vsk.line(grid.left, cell.position.y, grid.right, cell.position.y)


def draw_grid_border(vsk: vsketch.Vsketch, grid: Grid) -> None:
    vsk.rect(
        grid.top_left.x,
        grid.top_left.y,
        grid.width,
        grid.height,
    )
