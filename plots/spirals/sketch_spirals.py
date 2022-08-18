import vsketch
import random
import numpy as np
from typing import Sequence
from plotterdays.lib.constants.page_sizes import LETTER_SIZE_DIMENSIONS_METRIC
from shapely.geometry import LineString, Point
from shapely.affinity import scale
from plotterdays.lib.grid import Grid
from plotterdays.lib.spiral import Spiral


class SpiralsSketch(vsketch.SketchClass):
    # Sketch parameters:
    border = vsketch.Param(3.0)

    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("letter", landscape=True)
        vsk.scale("cm")

        # Setup
        num_rows = 10
        num_cols = 10
        grid_width = LETTER_SIZE_DIMENSIONS_METRIC[0] - (2 * self.border)
        grid_height = LETTER_SIZE_DIMENSIONS_METRIC[1] - (2 * self.border)

        grid = Grid(
            grid_width,
            grid_height,
            Point(self.border, self.border),
            num_cols,
            num_rows,
        )

        spiral_center = grid.center
        space_per_loop = 0.07

        # Draw the first spiral
        vsk.stroke(1)
        spiral_1 = Spiral(spiral_center, 0.0, 3.0, space_per_loop, 0, np.radians(1))
        spiral_1_points = spiral_1.get_points()
        self.draw_spiral(vsk, spiral_1_points, grid)

        # Draw the second spiral (& reflect it)
        vsk.stroke(2)
        spiral_2 = Spiral(spiral_center, 0.0, 3.0, space_per_loop, 0, np.radians(1))
        spiral_2_points = spiral_2.get_points()
        self.draw_spiral(vsk, spiral_2_points, grid, reflect=True)

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")

    def draw_spiral(
        self,
        vsk: vsketch.Vsketch,
        spiral_points: Sequence[Point],
        grid: Grid,
        reflect: bool = False,
    ) -> None:
        # This controls how spike-y / dense the noise is. The lower the value, the
        # more smooth the noise is (because we sample from a smaller range)
        spiral_noise_density = 500

        # TODO: Noise is a little too crazy right now.
        # Add some noise to the spiral
        x_coords = list(map(lambda v: v.x, spiral_points))
        y_coords = list(map(lambda v: v.y, spiral_points))

        noise_x_start = random.randrange(0, spiral_noise_density * 100)
        noise_x_space = vsk.noise(
            np.linspace(noise_x_start, spiral_noise_density, len(spiral_points)),
        )

        noise_y_start = random.randrange(10000, 10000 + (spiral_noise_density * 100))
        noise_y_space = vsk.noise(
            np.linspace(noise_y_start, noise_y_start + spiral_noise_density, len(spiral_points)),
        )

        coords_with_noise = []
        for i in range(len(x_coords)):
            noise_x = (noise_x_space[i] - 0.5) / 5
            noise_y = (noise_y_space[i] - 0.5) / 5
            noise_point = (x_coords[i] + noise_x, y_coords[i] + noise_y)
            coords_with_noise.append(noise_point)

        geom = LineString(coords_with_noise)
        if reflect:
            geom = scale(geom, xfact=-1.0, yfact=-1.0, origin=grid.center)

        vsk.geometry(geom)


if __name__ == "__main__":
    SpiralsSketch.display()
