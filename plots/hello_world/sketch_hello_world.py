import vsketch

class HelloWorldSketch(vsketch.SketchClass):
    # Sketch parameters:
    radius = vsketch.Param(2.0)

    def draw(self, vsk: vsketch.Vsketch) -> None:
        vsk.size("letter", landscape=True, center=False)
        vsk.penWidth("0.3mm")
        vsk.strokeWeight(3)
        # vsk.scale("cm")

        # implement your sketch here
        vsk.text("hello world!", 40.0, 100.0, size = 100.0, font='astrology', width = 800.0, align = "left")
        vsk.line(35.0, 160.0, 800.0, 160.0)
        vsk.line(800.0, 160.0, 35.0, 200.0)
        vsk.line(35.0, 200.0, 800.0, 200.0)

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")


if __name__ == "__main__":
    HelloWorldSketch.display()
