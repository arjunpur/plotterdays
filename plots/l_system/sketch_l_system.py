import numpy as np
from collections import deque
from lib.point_utils import add
from lib.vector import Vector, rotate_vector
import random
import vsketch
from typing import Tuple
from vsketch.shape import Point


class AnglePicker():
    def __init__(self, angle: float):
        self._angle = angle

    def angle(self) -> float:
        return self._angle


class RangeAnglePicker(AnglePicker):
    def __init__(self, angle_range: Tuple):
        self._range = angle_range

    def angle(self) -> float:
        return random.uniform(self._range[0], self._range[1])


class LSystemSketch(vsketch.SketchClass):
    # Sketch parameters:
    # radius = vsketch.Param(2.0)

    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("9in", "11in", landscape=True)
        vsk.scale("in")

        # system = LSystem(axiom="F", rules={"F": "FF+F-FF+FF"}, iterations=4)
        # draw_l_system(vsk, system, 60, 0.05, Point(0, 0))

        # Moving tree
        # system = LSystem(axiom="F+F", rules={"F": "FF[-FF][+F]"}, iterations=4)
        # draw_l_system(vsk, system, AnglePicker(30), 0.05, Point(0, 0), -90)

        # Stochastic moving tree
        # system = LSystem(axiom="F+F", rules={"F": "FF[-FFF][+F]"}, iterations=4)
        # draw_l_system(vsk, system, RangeAnglePicker((30, 45)), 0.05, Point(0, 0), -90)

        # system = LSystem(axiom="F", rules={"F": "FF[-FF][+FF]"}, iterations=4)
        # draw_l_system(vsk, system, RangeAnglePicker((0, 180)), 0.1, Point(0, 0), -90)

        system = LSystem(axiom="F", rules={"F": "FF[-FF][+FF]F"}, iterations=3)
        draw_l_system(vsk, system, AnglePicker(45), 0.2, Point(0, 0), -90)

        # implement your sketch here
        # vsk.circle(0, 0, self.radius, mode="radius")

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")


class LSystem():
    def __init__(self, axiom: str = "F", rules: object = {"F": "F+F-F-F+F"}, iterations: int = 2) -> None:
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


# Write a function that takes in an LSystem and draws it


def draw_l_system(
        vsk: vsketch.Vsketch,
        l_system: LSystem,
        angle_picker: AnglePicker,
        length: float,
        start: Point,
        starting_angle: float = 90,
):
    stack = deque()
    current_vector = Vector.from_length_and_angle(length, np.radians(starting_angle))
    current_position = start
    current_angle = 0
    for char in l_system.evaluate_production_rules():
        if char == "F":
            current_vector = rotate_vector(current_vector, current_angle, degrees=True)
            end_position = add(current_position, current_vector)

            vsk.line(current_position.x, current_position.y, end_position.x, end_position.y)

            current_position = end_position
            current_angle = 0
        elif char == "+":
            current_angle += angle_picker.angle()
        elif char == "-":
            current_angle -= angle_picker.angle()
        elif char == "[":
            stack.append((current_vector, current_position, current_angle))
        elif char == "]":
            current_vector, current_position, current_angle = stack.pop()


if __name__ == "__main__":
    LSystemSketch.display()
