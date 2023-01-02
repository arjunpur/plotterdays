from lib.curves import create_curve_with_light_bend, draw_cubic_bezier
from lib.grid import create_grid_with_padding, draw_grid
from lib.vector import Vector, draw_vector
import vsketch

class CurvesSketch(vsketch.SketchClass):
    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("9in", "11in", landscape=True, center=False)
        vsk.scale("in")

        grid = create_grid_with_padding(0.5, (0.50, 0.50)) 
        draw_grid(vsk, grid)

        rect_cell = grid.get_cell_at_index(0, 0)

        start = grid.get_cell_at_index(0, 0)
        end = grid.get_cell_at_index(6, 6)

        draw_vector(vsk, Vector.from_two_points(start.center, end.center).normalize(), start.center, draw_tip=True)

        curve1 = create_curve_with_light_bend(
            (start.center, end.center),
            (0.25, 0.8),
            start.width * 3,
            0.0,
            0.1,
            bend_clockwise = True,
        )
        draw_cubic_bezier(vsk, curve1)

        curve2 = create_curve_with_light_bend(
            (end.center, grid.get_cell_at_index(12, 12).center),
            (0.25, 0.8),
            start.width * 3,
            0.0,
            0.1,
            bend_clockwise = False,
        )
        draw_cubic_bezier(vsk, curve2)



