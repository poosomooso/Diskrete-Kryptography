import random
from elliptic_curve import e_curve, inverse_mod
# Pseudo-random number generator.
# Sequence looks like: f(x), f(f(x)), f(f(f(x))),...
# x - previous number
# a - constant of addition (wrote as parameter to be able to vary and compare search times)
# m - number to mod sequence by
def generate_next_elem(x, a, m):
	return (x**2 + a) % m

def calc_key(a, b, A, B, n):
	# key = (a - A) / (B - b) mod n
	return (a - A) * inverse_mod(B - b, n)

def run_pollard_rho_attack(curve, P, Q):
	tortoise_a = 1
	tortoise_b = 1
	hare_a = 1
	hare_b = 1

	# The algorithm should converge within around sqrt(n), where n is the order of the eliptic curve
	# Just in case, loop for n times. If it does not converge within this time, must try another generator
	print("Tortoise: (a, b) -> (x, y)  Hare: (A, B) -> (X, Y)")
	for j in range(curve.n):
		# Calculate aP + bQ for the tortoise and the hare random numbers
		tortoise_point = curve.add_points(curve.mult_point(P, tortoise_a), curve.mult_point(Q, tortoise_b))
		hare_point = curve.add_points(curve.mult_point(P, hare_a), curve.mult_point(Q, hare_b))

		print("Tortoise:", tortoise_point, "->", tortoise_point, " Hare:", hare_point, "->", hare_point)
		if(tortoise_point == hare_point):
			if(tortoise_a == hare_a):
				return -1 # No pairing found, key cannot be broken using this combination of generator/start points
			return calc_key(tortoise_a, tortoise_b, hare_a, tortoise_b, curve.n)
		tortoise_a = generate_next_elem(tortoise_a)
		tortoise_b = generate_next_elem(tortoise_b)
		hare_a = generate_next_elem(generate_next_elem(hare_a))
		hare_b = generate_next_elem(generate_next_elem(hare_b))

k = random.randint(1, curve.n)
run_pollard_rho_attack(e_curve, e_curve.g, e_curve.mult_point(e_curve.g, k))
