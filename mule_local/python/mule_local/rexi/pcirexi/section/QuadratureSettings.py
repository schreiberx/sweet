class QuadratureSettings:
    overall_quadrature_points: int = 1
    quadrature_method: str = "midpoint"
    distribute_quad_points_based_on_arc_length: bool = False

    def __init__(self, overall_quadrature_points, quadrature_method, distribute_quad_points_based_on_arc_length):
        self.overall_quadrature_points = overall_quadrature_points
        self.quadrature_method = quadrature_method
        self.distribute_quad_points_based_on_arc_length = distribute_quad_points_based_on_arc_length
