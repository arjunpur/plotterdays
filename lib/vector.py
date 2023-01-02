from __future__ import annotations
import numpy as np
import vsketch
import vpype
import latextools
from math import cos, sin
from shapely.geometry import Point, LineString, MultiLineString
from shapely.affinity import translate
from enum import Enum
from typing import Optional, Tuple
from io import StringIO

# TODO currently we proxy back and forth between numpy arrays
# I might be able to just remove this abstraction to simplify the code
# Or just use numpy arrays under the hood
class Vector:
    def __init__(self, x: float, y: float, scale: Tuple = (1.0, 1.0)):
        self.x = x * scale[0]
        self.y = y * scale[1]
        self.numpy = np.array([x, y])

    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.numpy))

    # https://stackoverflow.com/questions/44640479/type-annotation-for-classmethod-returning-instance
    @classmethod
    def from_numpy_array(cls, np_array: np.ndarray) -> Vector:
        assert(np_array.shape == (2,))
        return cls(np_array[0], np_array[1])

    @classmethod
    def from_two_points(cls, p1: Point, p2: Point) -> Vector:
        return Vector(p2.x - p1.x, p2.y - p1.y)

    def to_line_string(self, start_point: Point) -> LineString:
        end_point = Point(start_point.x + self.x, start_point.y + self.y)
        return LineString([start_point, end_point])

    def to_latex(self) -> str:
        # We need double braces because Python will confuse the brace with
        # string formatting.
        vector_str = r'''$
            \begin{{bmatrix}}
                {:0.2f} \\
                {:0.2f} \\
            \end{{bmatrix}}
        $'''.format(self.x, self.y).strip()
        return vector_str

    def get_perpendicular(self, clockwise = True) -> Vector:
        """
        Returns a vector perpendicular to this vector.
        Note: since the coordinate system's origin is in the top left, 
        the formula for calculating the perpendicular vector is inverted.
        """
        perp = self.numpy.copy()
        if clockwise:
            perp[0], perp[1] = -1 * self.y, self.x
        else:
            perp[0], perp[1] = self.y, -1 * self.x
        perp = perp / self.length
        return Vector(perp[0], perp[1])


    def normalize(self) -> Vector:
        if self.length == 0:
            raise Exception("could not compute norm of a 0 vector")
        return Vector.from_numpy_array(self.numpy / self.length)


    def dot_product(self, other: Vector) -> float:
        return np.dot(self.numpy, other.numpy)

    def scale(self, x, y) -> Vector:
        self.x *= x 
        self.y *= y 
        self.numpy = np.array([self.x, self.y])
        return self

    def __mul__(self, other: int | float) -> Vector:
        return Vector(self.x * other, self.y * other)

# TODO: Understand how vector rotation actually works mathematically
# TODO: Use shapely to perform the affine transformation. There's no
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
    rotated = np.dot(vec.numpy, rot_matrix)
    return Vector.from_numpy_array(rotated)


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

VECTOR_TIP_SCALING_FACTOR = -0.08
VECTOR_TIP_ANGLE_AWAY_FROM_MAIN = 20.0
def draw_vector(
    vsk: vsketch.Vsketch,
    vec: Vector,
    starting_point: Point,
    draw_tip: bool = True,
    draw_text: Optional[VectorTextOptions] = None,
):
    end_point = Point(starting_point.x + vec.x, starting_point.y + vec.y)
    vsk.line(starting_point.x, starting_point.y, end_point.x, end_point.y)
    
    # Create the tip of the vector by taking the original vector, rotating
    # in either direction + scaling it down.
    if draw_tip:
        one_side = rotate_vector(vec, VECTOR_TIP_ANGLE_AWAY_FROM_MAIN, degrees=True).normalize()
        other_side = rotate_vector(vec, -VECTOR_TIP_ANGLE_AWAY_FROM_MAIN, degrees=True).normalize()

        scaled_one_side = Vector(one_side.x * VECTOR_TIP_SCALING_FACTOR, one_side.y * VECTOR_TIP_SCALING_FACTOR)
        scaled_other_side = Vector(other_side.x * VECTOR_TIP_SCALING_FACTOR, other_side.y * VECTOR_TIP_SCALING_FACTOR)

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
        shapely_text.bounds
        
        # Plot the vector representation perpendicular to the vector. Adjust
        # the position too
        perp = vec.get_perpendicular() * 0.2
        if draw_text.orientation == VectorTextOrientation.LEFT_OF_VECTOR:
            perp = vec.get_perpendicular(clockwise = False) * 0.2

        mid_point = (starting_point.x + vec.x / 2, starting_point.y + vec.y / 2)
        text_location = Point(mid_point[0] + perp.x, mid_point[1] + perp.y)

        # Position the text to be at the tip of a scaled perpendicular
        # by finding the middle point of the text shape and translating the 
        # top left corner of shape to the tip of the vector and then adjusting
        # for the bounding box
        text_height = 0.5 * (shapely_text.bounds[3] - shapely_text.bounds[1])
        text_width = 0.5 * (shapely_text.bounds[2] - shapely_text.bounds[0])
        translated = translate(shapely_text, text_location.x - text_width, text_location.y - text_height)

        # TODO: TRANSLATE the text to the appropriate position based on 
        vsk.geometry(translated)


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
    

