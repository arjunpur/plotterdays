import numpy as np
import vsketch
from typing import Union
from shapely.geometry import Point
from plotterdays.lib.grid import Grid, Vector, draw_grid, draw_grid_border, rotate_vector, normalize
from plotterdays.lib.constants.page_sizes import LETTER_SIZE_DIMENSIONS_IMPERIAL


class VectorBasics(vsketch.SketchClass):
    # Sketch parameters:
    # radius = vsketch.Param(2.0)

    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("letter", landscape=True, center=False)
        vsk.scale("in")

        # Define the bounds of the sketch
        # AND Create the virtual grid
        #
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
        vsk.stroke(1)
        draw_grid_border(vsk, grid)

        draw_grid(vsk, grid)

        # Draw the vectors
        vec1 = Vector(1.0, 1.0)
        vec2 = Vector(2.0, 2.0)

        vsk.stroke(2)
        vec1 = grid.get_vector(2.0, 2.0)
        vec2 = grid.get_vector(-2.0, 2.0)
        starting_point = grid.get_cell_at_index(3, 3).center
        draw_vector_at(vsk, vec1, starting_point)
        draw_vector_at(vsk, vec2, starting_point)

        vsk.circle(starting_point.x, starting_point.y, 0.1)

        self.draw_star(vsk)

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


def draw_vector_at(vsk: vsketch.Vsketch, vec: Vector, point: Point, tip=True):
    end_point = Point(point.x + vec.x, point.y + vec.y)
    vsk.line(point.x, point.y, end_point.x, end_point.y)

    # Take the original vector, rotate it 20.0 degrees and -20.0, scale it's
    # normal by -0.1 to reverse and shorten it. This produces the vector's tip
    if tip:
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

if __name__ == "__main__":
    VectorBasics.display()
