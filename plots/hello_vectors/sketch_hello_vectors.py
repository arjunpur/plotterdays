from lib.grid import Grid, draw_grid, draw_grid_border
from lib.vector import Vector, draw_vector, VectorTextOptions, VectorTextOrientation
from lib.constants.page_sizes import LETTER_SIZE_DIMENSIONS_IMPERIAL

import vsketch
from shapely.geometry import Point

# TODO: I need to re-think how this entire thing is parametrized BytesIO
# scale. Currently the cell_size dictates the scale but the Vector
# abstraction doesn't really know anything about scale.
class VectorBasics(vsketch.SketchClass):
    def setup_grid(self) -> Grid:
        # NOTE: All units are in the imperial system (inches). This is because
        # the grid paper I'm using has well rounded units only in inches.
        padding = 0.5  # inches
        cell_size = (0.25, 0.25)
        # NOTE: we switch the access of width / height because we're doing landscape
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

        # Set up basic scale / font page_sizes
        font_size = grid.get_cell_at_index(0, 0).height

        # Draw the vectors
        vsk.stroke(2)
        cell_for_scale = grid.get_cell_at_index(0, 0)
        vec1 = Vector(2.0, 2.0, scale = (cell_for_scale.width, cell_for_scale.height))
        vec2 = Vector(-2.0, 2.0, scale = (cell_for_scale.width, cell_for_scale.height))
        starting_pos = (3, 3)
        starting_point = grid.get_cell_at_index(starting_pos[0], starting_pos[1]).center
        draw_vector(
            vsk,
            vec1,
            starting_point,
            draw_text=VectorTextOptions(VectorTextOrientation.RIGHT_OF_VECTOR, font_size),
        )

        draw_vector(
            vsk,
            vec2,
            starting_point,
            draw_text=VectorTextOptions(VectorTextOrientation.LEFT_OF_VECTOR, font_size),
        )

        vsk.circle(starting_point.x, starting_point.y, 0.1)

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")


# This is a different way of rendering LaTeX to SVGs, but unfortunately it 
# doesn't play well with svg
# def tex2svg(formula: str, fontsize: int = 30, dpi: int = 300) -> TextIOWrapper:
#     """Render TeX formula to SVG.
#     Args:
#         formula (str): TeX formula.
#         fontsize (int, optional): Font size.
#         dpi (int, optional): DPI.
#     Returns:
#         str: SVG render.
#     """

#     fig = plt.figure(figsize=(0.01, 0.01), frameon=False)
#     fig.text(0, 0, r'${}$'.format(formula), fontsize=fontsize)

#     output = BytesIO()
#     with open('test.svg', 'w+') as f:
#         fig.savefig(f, dpi=dpi, transparent=True, format='svg',
#                     bbox_inches='tight', pad_inches=0.0, frameon=False)
#         plt.close(fig)

#     output.seek(0)
#     return TextIOWrapper(output)

if __name__ == "__main__":
    VectorBasics.display()
