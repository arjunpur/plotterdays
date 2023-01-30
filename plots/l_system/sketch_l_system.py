from lib.point_utils import add
from lib.vector import Vector, rotate_vector
import vsketch
from vsketch.shape import Point

class LSystemSketch(vsketch.SketchClass):
    # Sketch parameters:
    # radius = vsketch.Param(2.0)

    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("9in", "11in", landscape=True)
        vsk.scale("in")

        system = LSystem(axiom="F", rules={"F": "FF+F-FF+FF"}, iterations=4)
        draw_l_system(vsk, system, 60, 0.05, Point(0, 0))

        # implement your sketch here
        # vsk.circle(0, 0, self.radius, mode="radius")

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")


class LSystem():
    def __init__(self, axiom="F", rules={"F": "F+F-F-F+F"}, iterations=2):
        self.axiom = axiom
        self.rules = rules
        self.iterations = iterations
        self.output = self.evaluate_production_rules()


    def evaluate_production_rules(self) -> str:
        output = self.axiom
        for _ in range(self.iterations):
            new_output = ""
            for char in output:
                if char in self.rules:
                    new_output += self.rules[char]
                else:
                    new_output += char
            output = new_output
        return output


def draw_l_system(vsk: vsketch.Vsketch, l_system: LSystem, angle: float, length: float, start: Point):
    current_vector = Vector.from_two_points(start, Point(start.x + length, start.y + length))
    current_position = start
    current_angle = 0 
    for char in l_system.evaluate_production_rules():
        if char == "F":
            current_vector = rotate_vector(current_vector, current_angle, degrees = True)
            end_position = add(current_position, current_vector)

            vsk.line(current_position.x, current_position.y, end_position.x, end_position.y)

            current_position = end_position 
            current_angle = 0
        elif char == "+":
           current_angle += angle
        elif char == "-":
            current_angle -= angle

if __name__ == "__main__":
    LSystemSketch.display()
