from typing import Tuple
from lib.constants.page_sizes import LETTER_SIZE_DIMENSIONS_IMPERIAL, LETTER_SIZE_DIMENSIONS_METRIC
import vsketch
from shapely.geometry import Point

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

    @property
    def dimensions(self) -> Tuple:
        return (self.width, self.height) 


class Grid:
    """
    A Cartesian grid that is agnostic of units. Callers are expected to maintain consistency of
    units when using the grid.

    `top_left` refers to the Grid's top-left coordinates.
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
        self.cell_width = self.width / self.num_cols
        self.cell_height = self.height / self.num_rows

    def get_cell_at_index(self, row: int, col: int) -> Cell:
        width = self.cell_width
        height = self.cell_height
        cell_x = self.top_left.x + (col * width)
        cell_y = self.top_left.y + (row * height)
        return Cell(Point(cell_x, cell_y), width, height)

    @property
    def cell_dimensions(self) -> Tuple:
        return self.get_cell_at_index(0, 0).dimensions

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


def create_grid_with_padding(padding: float, cell_size: Tuple, use_imperial: bool = True, landscape: bool = True) -> Grid:
    """
    Returns a Grid object who's coordinates have been adjusted by `padding`, with cell_size
    the size of `cell_size`. 
    The units of measurement (imperial / metric) and the orientation of the page (landscape, potrait) 
    should be specified.
    """
    if use_imperial:
        dimensions = LETTER_SIZE_DIMENSIONS_IMPERIAL
    else:
        dimensions = LETTER_SIZE_DIMENSIONS_METRIC
    # NOTE: we switch the access of width / height because we're doing landscape
    if landscape: 
        grid_width = dimensions[1] - (2 * padding)
        grid_height = dimensions[0] - (2 * padding)
    else: 
        grid_width = dimensions[0] - (2 * padding)
        grid_height = dimensions[1] - (2 * padding)
    num_cols = int(round(grid_width / cell_size[1]))
    num_rows = int(round(grid_height / cell_size[0]))
    grid = Grid(
        grid_width,
        grid_height,
        Point(padding, padding),
        num_cols,
        num_rows,
    )
    return grid

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
