from plotterdays.lib.grid import Grid, Vector, draw_grid, draw_grid_border, rotate_vector, normalize
from plotterdays.lib.constants.page_sizes import LETTER_SIZE_DIMENSIONS_IMPERIAL

import numpy as np
import vsketch
from typing import Optional, Tuple
from enum import Enum
from shapely.geometry import Point
import vpype
from io import StringIO
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


class VectorBasics(vsketch.SketchClass):
    # Sketch parameters:
    # radius = vsketch.Param(2.0)

    def setup_grid(self) -> Grid:
        # NOTE: All units are in the imperial system (inches). This is because
        # the grid paper I'm using has well rounded units only in inches.
        padding = 0.5  # inches
        cell_size = (0.25, 0.25)
        # NOTE: we switch the order of width / height because we're doing landscape
        grid_width = LETTER_SIZE_DIMENSIONS_IMPERIAL[1] - (2 * padding)
        grid_height = LETTER_SIZE_DIMENSIONS_IMPERIAL[0] - (2 * padding)
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

    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("letter", landscape=True, center=False)
        vsk.scale("in")

        # Define the bounds of the sketch
        # AND Create the virtual grid
        grid = self.setup_grid()
        vsk.stroke(1)
        draw_grid_border(vsk, grid)
        draw_grid(vsk, grid)

        # Draw the vectors
        vsk.stroke(2)
        vec1 = grid.get_vector(2.0, 2.0)
        vec2 = grid.get_vector(-2.0, 2.0)
        starting_pos = (3, 3)
        draw_vector_at(
            vsk,
            grid,
            starting_pos,
            vec1,
            draw_tip=True,
            draw_text=VectorTextOrientation.RIGHT_OF_VECTOR,
        )
        draw_vector_at(vsk, grid, starting_pos, vec2)

        point = grid.get_cell_at_index(starting_pos[0], starting_pos[1]).center
        vsk.circle(point.x, point.y, 0.1)

        # self.draw_star(vsk)

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")


    def draw_star(self, vsk: vsketch.Vsketch) -> None:
        angles = np.linspace(0, 2 * np.pi, 5, endpoint=False)
        idx = [0, 2, 4, 1, 3, 0]
        x = np.cos(angles[idx] - np.pi / 2)
        y = np.sin(angles[idx] - np.pi / 2)

        with vsk.pushMatrix():
            for i in range(5):
                with vsk.pushMatrix():
                    vsk.scale(0.8**i)
                    vsk.polygon(x, y)

                vsk.translate(2, 0)

        vsk.translate(0, 4)

        for i in range(5):
            with vsk.pushMatrix():
                vsk.rotate(i * 4, degrees=True)
                vsk.polygon(x, y)

            vsk.translate(2, 0)


class VectorTextOrientation(Enum):
    RIGHT_OF_VECTOR = 1
    LEFT_OF_VECTOR = 2
    START_OF_VECTOR = 3
    END_OF_VECTOR = 4


def draw_vector_at(
    vsk: vsketch.Vsketch,
    grid: Grid,
    starting_pos: Tuple[int, int],
    vec: Vector,
    draw_tip: bool = True,
    draw_text: Optional[VectorTextOrientation] = None,
):
    starting_cell = grid.get_cell_at_index(starting_pos[0], starting_pos[1])
    point = starting_cell.center
    end_point = Point(point.x + vec.x, point.y + vec.y)
    vsk.line(point.x, point.y, end_point.x, end_point.y)

    # Take the original vector, rotate it 20.0 degrees and -20.0, scale it's
    # normal by -0.1 to reverse and shorten it. This produces the vector's tip
    if draw_tip:
        one_side = normalize(rotate_vector(vec, 20.0, degrees=True))
        other_side = normalize(rotate_vector(vec, -20.0, degrees=True))

        scaled_one_side = Vector(one_side.x * -0.1, one_side.y * -0.1)
        scaled_other_side = Vector(other_side.x * -0.1, other_side.y * -0.1)

        vsk.line(
            end_point.x,
            end_point.y,
            end_point.x + scaled_one_side.x,
            end_point.y + scaled_one_side.y
        )

        vsk.line(
            end_point.x,
            end_point.y,
            end_point.x + scaled_other_side.x,
            end_point.y + scaled_other_side.y
        )
    if draw_text is not None:
        pass
        # Part 1
        # Generate the Latex SVG

        vector_repr_svg = tex2svg(r'E=mc^2')
        vpype_text, _, _ = vpype.read_svg(vector_repr_svg, 0.001)

        # Part 2
        # Position it / scale it to the grid / vector / units
        # Convert to MultiLineString using vpype
        # Sketch the shape using vsketch
        # The text should only be the size of two cells

        scale_factor = (starting_cell.height * 2.0) / vpype_text.height()
        vpype_text.scale(scale_factor)
        vpype_text.translate(starting_cell.left, starting_cell.top)
        shapely_text = vpype_text.as_mls()
        vsk.geometry(shapely_text)


def tex2svg(formula: str, fontsize: int = 20, dpi: int = 300) -> StringIO:
    """Render TeX formula to SVG.
    Args:
        formula (str): TeX formula.
        fontsize (int, optional): Font size.
        dpi (int, optional): DPI.
    Returns:
        str: SVG render.
    """

    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, r'${}$'.format(formula), fontsize=fontsize)

    output = StringIO()
    fig.savefig(output, dpi=dpi, transparent=True, format='svg',
                bbox_inches='tight', pad_inches=0.0, frameon=False)
    plt.close(fig)

    output.seek(0)
    return output


if __name__ == "__main__":
    VectorBasics.display()
