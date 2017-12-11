from elliptic_curve import EllipticCurve

import math

def bsgs(curve):
	"""
	Idea is to find x such that the public key is x * base on the 
	given elliptic curve.

	Since all integers can be represented as x = am + b, we break down the search space
	by searching for a and b, which are both less than sqrt(n)
	"""
	base = curve.g
	pub_key = curve.get_pub_key()

	m = int(math.sqrt(curve.n)) # sqrt of the order of the group

	baby_step_table = {}

	for b in range(m):
		baby = curve.mult_point(base, b)
		baby_step_table[baby] = b

	step = curve.mult_point(base, m) # m * base
	step = curve.neg_point(step)

	required_baby = pub_key # pub-key - a*m*base

	for a in range(m):
		baby_size = baby_step_table.get(required_baby)

		if baby_size:
			return baby_size + (a * m)

		required_baby = curve.add_points(required_baby, step) # subtract another m*base

def test_bsgs():
	curve = EllipticCurve(10177, 1, -1, 10331, (1, 1))
	priv_key = bsgs(curve)

	print('calculated key:', priv_key, 'actual key:', curve._d, 'is equal:', priv_key == curve._d)

if __name__ == '__main__':
	test_bsgs()