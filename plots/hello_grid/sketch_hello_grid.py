import vsketch


class HelloGridSketch(vsketch.SketchClass):
    # Sketch parameters:
    # radius = vsketch.Param(2.0)
    grid_rows = vsketch.Param(10.0)
    grid_cols = vsketch.Param(10.0)
    grid_width = vsketch.Param(10.0)
    grid_height = vsketch.Param(10.0)

    def draw(self, vsk: vsketch.Vsketch) -> None:
        # Center ensures the top left coordinates of the page are (0, 0) with 
        # coordinates increasing as you go down
        vsk.size("letter", landscape=True, center=False)
        vsk.scale("cm")

        page_size = vsk.document.page_size
        print(page_size)
        print(vsk.document.metadata)

        # Draw bounding rectangle
        # Fill lines with grid resolution
        # implement your sketch here
        # vsk.circle(0, 0, self.radius, mode="radius")

    def finalize(self, vsk: vsketch.Vsketch) -> None:
        vsk.vpype("linemerge linesimplify reloop linesort")


if __name__ == "__main__":
    HelloGridSketch.display()
