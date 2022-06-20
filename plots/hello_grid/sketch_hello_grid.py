import vsketch
import numpy as np
from shapely.geometry import LineString
from typing import Tuple


class HelloGridSketch(vsketch.SketchClass):
    # Sketch parameters:
    page_dimensions = (28.0, 21.5)
    num_rows = vsketch.Param(10.0)
    num_cols = vsketch.Param(10.0)
    border = vsketch.Param(3.0)

    def draw(self, vsk: vsketch.Vsketch) -> None:
        # Center ensures the top left coordinates of the page are (0, 0) with
        # coordinates increasing as you go down
        # Dimensions: 28cm by 21.5cm
        vsk.size("letter", landscape=True, center=False)
        vsk.scale("cm")
        vsk.detail("0.1mm")

        page_size = vsk.document.page_size
        print(page_size)
        print(vsk.document.metadata)

        # TODO: Parametrize this grid
        # self.draw_grid(vsk)

        grid_dimensions = (
            self.page_dimensions[0] - (2 * self.border),
            self.page_dimensions[1] - (2 * self.border),
        )
        col_increment = grid_dimensions[0] / self.num_cols
        row_increment = grid_dimensions[1] / self.num_rows

        spiral_center = (self.border + col_increment * 5.0, self.border + row_increment * 5.0)
        self.draw_spiral(vsk, spiral_center, 0.0, 0.1, 0.0, 10 * np.pi, np.radians(1))

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")

    def draw_grid(self, vsk: vsketch.Vsketch) -> None:
        grid_dimensions = (
            self.page_dimensions[0] - (2 * self.border),
            self.page_dimensions[1] - (2 * self.border),
        )

        vsk.rect(
            self.border,
            self.border,
            grid_dimensions[0],
            grid_dimensions[1],
        )

        col_increment = grid_dimensions[0] / self.num_cols
        row_increment = grid_dimensions[1] / self.num_rows

        for row in range(1, int(self.num_rows)):
            # Draw the column lines (vertical)
            for col in range(1, int(self.num_cols)):
                column_x = self.border + (col * col_increment)
                column_y_start = self.border
                column_y_end = self.border + grid_dimensions[1]
                vsk.line(column_x, column_y_start, column_x, column_y_end)
            row_y = self.border + (row * row_increment)
            row_x_start = self.border
            row_x_end = self.border + grid_dimensions[0]
            vsk.line(row_x_start, row_y, row_x_end, row_y)

    # TODO: Change the interface of this to take in start_radius and end_radius
    # instead of start_theta and end_theta. It's more ergonomical to use distances
    # than angles to express the size of the spiral.
    def draw_spiral(
        self,
        vsk: vsketch.Vsketch,
        center: Tuple[float, float],
        start_radius: float,
        space_per_loop: float,
        start_theta: float,
        end_theta: float,
        theta_step: float,
    ) -> None:

        # This uses the basic formula for determining the radius to use
        # when drawing a point on a spiral (in polar coordinates)
        # The theta parameter will keep changing. `b` represents the space_per_loop
        # between each loop, and `a` represents where you're starting
        #
        # https://en.wikipedia.org/wiki/Archimedean_spiral
        def calculate_x_y(a: float, b: float, theta: float):
            r = a + (b * theta)
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            return (center[0] + x, center[1] + y)

        a = start_radius
        b = space_per_loop
        temp_theta = start_theta
        start_point = calculate_x_y(a, b, temp_theta)
        spiral_points = [start_point]

        while temp_theta < end_theta:
            temp_theta += theta_step
            new_point = calculate_x_y(a, b, temp_theta)
            spiral_points.append(new_point)

        geom = LineString(spiral_points)
        vsk.geometry(geom)
        pass


if __name__ == "__main__":
    HelloGridSketch.display()
