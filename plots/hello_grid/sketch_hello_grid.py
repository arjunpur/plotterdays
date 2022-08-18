import vsketch
from shapely.geometry import Point
from plotterdays.lib.grid import Grid, draw_grid


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

        width = self.page_dimensions[0] - (2 * self.border)
        height = self.page_dimensions[1] - (2 * self.border)
        grid = Grid(
            width,
            height,
            Point(self.border, self.border),
            int(self.num_cols),
            int(self.num_rows),
        )
        draw_grid(vsk, grid)

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")


if __name__ == "__main__":
    HelloGridSketch.display()
