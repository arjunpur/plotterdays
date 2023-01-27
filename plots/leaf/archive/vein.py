
### Archive
class Vein():

    def __init__(
        self,
        vein_start: Point,
        vein_trajectory: Vector,
        vein_boundaries: List[bezier.Curve],
        level: int = 0,
        vsk: Optional[vsketch.Vsketch] = None,
    ):
        self.level = level
        if vsk:
            self.vsk = vsk

        for boundary in vein_boundaries:
            computed_vein_end = self._compute_vein_end(vein_start, vein_trajectory, boundary)
            if computed_vein_end:
                self.vein_end = computed_vein_end
                break
        if not self.vein_end:
            raise Exception("could not compute vein end - invalid boundaries")

        vein_line = LineString([vein_start, self.vein_end])
        self.vein_boundary = vein_boundaries
        
        # TODO: This might not work if the orientation of the vein is reverse
        trajectory = Vector.from_two_points(vein_start, self.vein_end).normalize()
        left_or_right = trajectory.x > 0
        self.curve = create_curve_with_light_bend(
            (vein_start, self.vein_end),
            (0.25, 0.75),
            vein_line.length * 0.2,
            vein_line.length * 0.1,
            vein_line.length * 0.1,
            bend_clockwise = left_or_right,
        )
        
    def _compute_vein_end(self, vein_start: Point, vein_trajectory: Vector, side_curve: bezier.Curve) -> Optional[Point]:
        try:
            scaled_vein_trajectory = vein_trajectory * sys.maxsize
            bounding_point_end = add(vein_start, scaled_vein_trajectory)
            nodes = np.asarray([[vein_start.x, bounding_point_end.x], [vein_start.y, bounding_point_end.y]])
            line_segment = bezier.Curve(
                nodes,
                degree = 1
            )
            intersection_s_value_with_vein_boundary = side_curve.intersect(line_segment)[0][0]
            intersection_point = side_curve.evaluate(intersection_s_value_with_vein_boundary)
            return Point(intersection_point)
        except:
            return None
