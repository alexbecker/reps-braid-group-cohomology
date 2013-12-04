# CyclePolynomials module
# implements the CyclePoly class and related methods	

import Perm
from numpy import *
from numpy.linalg import lstsq

# zero cutoff for coefficients of a CyclePoly
epsilon = .0005

# returns all monomials of degree d in r cycle variables
# a monomial is represented as a list of integers, with the ith entry being the degree of the (i+1)st cycle variable
# this function is shit, but it's not a bottleneck so w/e
def monomials(d, r):
	if d == 0:
		return [[0 for i in range(r)]]

	result = []
	for monomial in monomials(d - 1, r):
		for i in range(r):
			result.append(monomial[:i] + [monomial[i] + 1] + monomial[i + 1:])

	result.sort()

	# nub result
	i = 0
	while i < len(result) - 1:
		if result[i] == result[i+1]:
			result.pop(i)
		else:
			i += 1

	return result

# assigns a weight to a monomial which is higher for larger cycle variables
# used to heuristically trim the search space of polynomials for interpolation
def weight(monomial):
	multiplier, result = 1, 0
	for degree in monomial:
		result += multiplier * degree
		multiplier += 3.3 / len(monomial)

	return result 

def evalMonomial(monomial, cycleType):
	result, i = 1, 1
	for degree in monomial:
		result *= cycleType.count(i) ** degree
		i += 1

	return result

# defines cycle polynomials as a list of monomials and a list of coefficients
class CyclePoly:
	def __init__(self, terms, coefficients):
		self.terms = terms
		self.coefficients = coefficients

	def __str__(self):
		result = ''
		for i in range(len(self.terms)):
			coeff = self.coefficients[i]
			if abs(coeff) > epsilon:
				result += str(coeff)
				j = 1
				for degree in self.terms[i]:
					if degree != 0:
						result += 'X_{}^{}'.format(j, degree)
					j += 1
				result += ' + '
		
		# return all but the last ' + '
		return result[:len(result) - 3]

	# return the monomials with substantially nonzero coefficients
	def nonzeroTerms(self):
		result = []
		for i in range(len(self.terms)):
			if abs(self.coefficients[i]) > epsilon:
				result.append(self.terms[i])

		return result

# given a list of values of a function f on the cycle types of S_n (sorted lexigraphically)
# returns the best degree d cycle polynomial interpolation of f in r cycle variables
def interpolateCyclePoly(n, functionVals, d, r):
	points = Perm.partitions(1, n)

	# all monomials of degree at most d
	terms = []
	for i in range(d+1):
		terms.extend(monomials(i, r))

	# let m_i be the monomials of degree at most d in r cycle variables, sorted lexigraphically
	# let c_i be the cycle types of S_n, sorted lexigraphically
	# A * x returns the value of sum_i x_im_i(c_j)
	A = empty( (len(points), len(terms)) )
	for i in range(len(points)):
		for j in range(len(terms)):
			A[i,j] = evalMonomial(terms[j], points[i])

	coefficients, lstsqError = lstsq(A, functionVals)[:2]

	print('Error: {}'.format(lstsqError))

	return CyclePoly(terms, coefficients)

# better version of interpolateCyclePolynomial for functions which are defined on all S_i
# given the list of values of a function f on the cycles types of S_x through S_n (sorted by S_i, then lexigraphically)
# returns the best cycle polynomials interpolation (using the given terms) of f
# works by minimizing the error across S_x through S_n
def multiInterpolate(x, n, functionVals, terms):
	points = []
	for k in range(x,n+1):
		points.extend(Perm.partitions(1, k))

	# like A from interpolateCyclePolynomials, but A * x returns the concatenation of all values sum_i x_im_i(c_j)
	A = empty( (len(points), len(terms)) )
	for i in range(len(points)):
		for j in range(len(terms)):
			A[i,j] = evalMonomial(terms[j], points[i])

	coefficients, lstsqError = lstsq(A, functionVals)[:2]

	print('Error: {}'.format(lstsqError))

	return CyclePoly(terms, coefficients)

