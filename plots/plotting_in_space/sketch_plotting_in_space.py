from plotterdays.lib.grid import Grid, Vector, draw_grid, draw_grid_border, rotate_vector, normalize
from plotterdays.lib.constants.page_sizes import LETTER_SIZE_DIMENSIONS_IMPERIAL

import latextools
import vsketch
from typing import Optional, Tuple
from enum import Enum
from shapely.geometry import Point, MultiLineString
import vpype
from io import StringIO

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
        vec1 = grid.get_vector(2.0, 2.0)
        vec2 = grid.get_vector(-2.0, 2.0)
        starting_pos = (3, 3)
        starting_point = grid.get_cell_at_index(starting_pos[0], starting_pos[1]).center
        draw_vector_at(
            vsk,
            vec1,
            starting_point,
            draw_text=VectorTextOptions(VectorTextOrientation.RIGHT_OF_VECTOR, font_size),
        )
        draw_vector_at(vsk, vec2, starting_point)

        point = grid.get_cell_at_index(starting_pos[0], starting_pos[1]).center
        vsk.circle(point.x, point.y, 0.1)

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")

class VectorTextOrientation(Enum):
    RIGHT_OF_VECTOR = 1
    LEFT_OF_VECTOR = 2
    START_OF_VECTOR = 3
    END_OF_VECTOR = 4

# TODO: Is this pythonic?
class VectorTextOptions:
    def __init__(self, orientation: VectorTextOrientation, font_size: float):
        self.orientation = orientation
        self.font_size = font_size


def draw_vector_at(
    vsk: vsketch.Vsketch,
    vec: Vector,
    starting_point: Point,
    draw_tip: bool = True,
    draw_text: Optional[VectorTextOptions] = None,
):
    end_point = Point(starting_point.x + vec.x, starting_point.y + vec.y)
    vsk.line(starting_point.x, starting_point.y, end_point.x, end_point.y)

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
        # Part 1 - Generate the Latex SVG
        
        shapely_text = latex2shapely(vec.to_latex(), draw_text.font_size)

        # TODO: TRANSLATE the text to the appropriate position based on 
        vsk.geometry(shapely_text)


#TODO Abstract this function further so that it returns a Shapely / vpype object
def latex2shapely(formula: str, font_size: float = 1.0, debug: bool = False) -> MultiLineString:
    pdf = latextools.render_snippet(
        formula,
        commands=[latextools.cmd.all_math],
        pad = '0.001in' # Pad the latex SVG to avoid it getting truncated
    )
    svg_eq = pdf.as_svg().content
    if debug:
        with open('latex2svg.svg', 'w+') as f:
            f.write(svg_eq)

    svg_eq_buffer = StringIO(svg_eq)
    vpype_text, _, _ = vpype.read_svg(
        svg_eq_buffer,
        0.001, # quantization - how granular the lines are for curves
        simplify = True,
    )

    scale_factor = font_size / vpype_text.height()
    vpype_text.scale(scale_factor)

    # Align the SVG with the top left corner of the page 
    bounding_box = vpype_text.bounds()
    if bounding_box is None:
        raise Exception("unexpectedly encountered None bounding_box for LaTex SVG")
    to_translate_by = (bounding_box[0], bounding_box[1])
    vpype_text.translate(-1 * to_translate_by[0], -1 * to_translate_by[1])

    return vpype_text.as_mls()
    

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
