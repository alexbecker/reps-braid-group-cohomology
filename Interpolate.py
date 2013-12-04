# Interpolate module
# interpolates the cycle polynomial for the character from the bonus problem from Math 267, HW 6
# Usage: python Interpolate.py i filename

import CyclePolynomials as cp
import pickle
import numpy
from json import loads
from Bonus import *

# NOTE: heuristic, not rigorous
# returns the most likely monomials of degree at most d in r cycle variables to appear in the character's cycle polynomial
def likelyMonomials(d, r):
	allMonomials = []
	for i in range(1,d+1):
		allMonomials.extend(cp.monomials(i, r))

	return [monomial for monomial in allMonomials if cp.weight(monomial) <= d]

# interpolates the best cycle polynomial for the ith desired character
# works by minimizing the error across V_x through V_n where x is as small as possible s.t. its ith exterior product is nonzero
# obtains the degree d and number of variables r from stdin, or alternatively a list of monomials
# retrieves the character values from filename
def cyclePolynomial(i, filename):
	# compute the smallest possible x such that the ith exterior power of x is nonzero
	x = 2
	while len(V(x)) < i:
		x += 1

	# retrieve character values
	inputFile = open(filename, 'rb')
	inputFile.seek(0, 2)
	fileLength = inputFile.tell()
	inputFile.seek(0)
	characterVals = []
	while inputFile.tell() < fileLength:
		characterVals.append(pickle.load(inputFile))

	inputFile.close()

	# sort charcterVals just in case
	characterVals.sort()

	# filter character values to get characters for correct i and remove characters for V_k, k < x
	characterVals = filter(lambda y: y[0][0] >= x and y[0][1] == i, characterVals)

	# adjust x to be the smallest k such that the character of V_k is included in characterVals
	x = characterVals[0][0][0]	
	# set n to be the largest
	n = characterVals[len(characterVals) - 1][0][0]

	# remove the tuples from characterVals and flatten
	characterVals = map(lambda y: y[1], characterVals)
	characterVals = sum(characterVals, [])

	# get d and r (or all terms) from stdin and interpolate until user is satisfied
	# user can opt to refine the result by restricting to terms with substantially nonzero coefficient
	satisfied = 'n'
	while satisfied != 'y':
		# ask user whether to supply terms
		supplyTerms = raw_input('Enter terms manually? (y/n): ')
		if supplyTerms == 'y':
			terms = loads(raw_input('Enter terms: '))
		else:
			d = int(raw_input('Enter d: '))
			r = int(raw_input('Enter r: '))
			terms = likelyMonomials(d, r)

		# find solution
		soln = cp.multiInterpolate(x, n, characterVals, terms)
		print(soln)

		# refine soln if desired
		refine = raw_input('Refine? (y/n): ')
		while refine == 'y':
			terms = soln.nonzeroTerms()
			soln = cp.multiInterpolate(x, n, characterVals, terms)
			print(soln)
			refine = raw_input('Refine? (y/n): ')

		satisfied = raw_input('Satisfied? (y/n): ')

	return soln

# MAIN: runs cyclePolynomial(n, i)
if __name__ == '__main__':
	cyclePolynomial(int(argv[1]), argv[2])

