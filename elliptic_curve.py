""" Modeled after https://github.com/andreacorbellini/ecc/blob/master/logs/common.py"""

import random

def inverse_mod(n, p):
    """
    Returns a n_inv such that (n * n_inv) === 1 modp
    n must be nonzero and p must be prime
    solved using the euclidean algorithm
    """
    if n == 0:
        raise ZeroDivisionError()
    if n < 0:
        # wrap around
        return p - inverse_mod(-n, p)

    curr_s, prev_s = 0, 1
    curr_r, prev_r = p, n

    while curr_r != 1:
        quotient = prev_r // curr_r
        prev_r, curr_r = curr_r, (prev_r - quotient * curr_r)
        prev_s, curr_s = curr_s, (prev_s - quotient * curr_s)

    gcd, x = curr_r, curr_s

    assert gcd == 1
    assert (n * x) % p == 1

    return x % p

class EllipticCurve:
    """ Equation of the form y^2=x^3+ax+b """
    """ A point being None means that it is the point at infinity (e.g 0) """
    def __init__(self, field, a, b, n, g):
        """
        field : size of the finite field
        a & b : parameters on the curve
        n : order of the cyclic subgroup (pregenerated using schoof's algorithm)
        g : the base point
        """
        self.field = field
        self.a = a
        self.b = b
        self.n = n
        self.g = g
        self._d = random.randrange(1, self.n) # private key, randomly generated

        # making sure the parameters are valid
        assert pow(2, field - 1, field) == 1
        assert (4 * a * a * a + 27 * b * b) % field != 0
        assert self.is_on_curve(g)
        assert self.mult_point(g, n) is None

    def get_pub_key(self):
        return self.mult_point(self.g, self._d)

    def gen_shared(self, other_point):
        self._shared = self.mult_point(self.g, self._d, starting_point=other_point)

    def is_on_curve(self, point):
        if point is None:
            return True

        x, y = point
        rhs = (x * x * x) + (self.a * x) + self.b
        lhs = y * y
        return (lhs - rhs) % self.field == 0

    def add_points(self, point1, point2):
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
            # same as finding the slope normally, but we need to use modular
            # inverse instead of multiplicative inverse
            m = (y1 - y2) * inverse_mod(x1 - x2, self.field)

        x3 = m * m - x1 - x2
        y3 = y1 + m * (x3 - x1)
        result = (x3 % self.field,
                  -y3 % self.field)

        assert self.is_on_curve(result)
        return result

    def neg_point(self, point):
        """ Returns the negative point, i.e. x, -y """
        if point is None:
            return None

        x, y = point
        result = x, -y % self.field

        assert self.is_on_curve(result)
        return result

    def mult_point(self, point, n, starting_point=None):
        """ returns n * point i.e. point added to itself n times """
        
        if n % self.n == 0 or point is None:
            return None

        if n < 0: # go the opposite way
            return self.neg_point(self.mult_point(point, -n))

        result = starting_point
        addend = point # addend is a fancy name for an addition operand

        while n:
            if n & 1:
                result = self.add_points(result, addend)

            # double the addend
            addend = self.add_points(addend, addend)

            # halve the number of points to multiply by
            n >>= 1

        return result

e_curve = EllipticCurve(
    field=10177,
    a=1,
    b=-1,
    g=(1, 1),
    n=10331,
)
