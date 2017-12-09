""" Modeled after https://github.com/andreacorbellini/ecc/blob/master/logs/common.py"""

class EllipticCurve:
    """ Equation of the form y^2=x^3+ax+b """
    """ A point being None means that it is the point at infinity (e.g 0) """
    def __init__(self):
        self.field = field
        self.a = a
        self.b = b

        # making sure the parameters are good
        assert pow(2, field - 1, field) == 1
        assert (4 * a * a * a + 27 * b * b) % field != 0

    def add_pub_keys(g, n):
        self.g = g
        self.n = n

        # make sure keys are valid
        assert self.is_on_curve(g)
        assert self.mult(n, g) is None

    def is_on_curve(self, point):
        if point is None:
            return True

        x, y = point
        rhs = (x * x * x) + (self.a * x) + self.b
        lhs = y * y
        return (lhs - rhs) % self.field == 0

    def add_points(point1, point2):
        assert self.is_on_curve(point1)
        assert self.is_on_curve(point2)
        
        # If a point is 0, return the other point
        if point1 is None:
            return point2
        if point2 is None:
            return point1

        x1, y1 = point1
        x2, y2 = point2

        if x1 == x2 and y1 != y2: 
            # these points are opposite; the line is vertical and we return 0
            return None

        # these next ones require us to find a slope, m
        if x1 == x2:
            # tangent line
            m = (3 * x1 * x1 + self.a) * inverse_mod(2 * y1, self.field)
        else:
            # line through both points
            m = (y1 - y2) * inverse_mod(x1 - x2, self.field)

        x3 = m * m - x1 - x2
        y3 = y1 + m * (x3 - x1)
        result = (x3 % self.field,
                  -y3 % self.field)

        assert self.is_on_curve(result)
        return result