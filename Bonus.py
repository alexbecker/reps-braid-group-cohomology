# Bonus module
# contains code specific for the Bonus problem from Math 267, HW 6
# USAGE: python Bonus.py m n i filename

import ExteriorAlg as ea
import Subsets as ss
import Perm
import pickle
from sys import argv
from scipy.sparse import *
from scipy.sparse.linalg import spsolve
from multiprocessing import Pool
from parmap import parmap

# helper function for various functions
def sgn(m):
	if m % 2 == 0:
		return 1
	return -1

# returns (i,j) if i<j and (j,i) if j<i
# NOTE: this is NOT being viewed as its first exterior power
def makePair(i, j):
	if i < j:
		return (i, j)
	else:
		return (j, i)

# returns a basis for the vector space V_n
# NOTE: indexed from 0 rather than 1
def V(n):
	return [(i, j) for i in range(n) for j in range(n) if i < j]

# applies a permutation to a basis vector of V_n
def applyPermV(perm, elt):
	i, j = elt
	return makePair(perm[i], perm[j])

# defines the vector R_{i,j,k}
def R(j, k, l):
	return ea.Element(list(map(ea.standardForm, [(1, [makePair(j, k), makePair(k, l)]), 
												 (1, [makePair(k, l), makePair(l, j)]),
												 (1, [makePair(l, j), makePair(j, k)])])))

# applies a permutation to an element of the kth exterior power of V_n
def applyPermExtV(perm, elt):
	permutedElt = []
	for coeff, basisVector in elt.coeffVectorPairs:
		permutedElt.append((coeff, list(map(lambda x: applyPermV(perm, x), basisVector))))

	result = ea.Element(permutedElt)
	result.standardForm()

	return result

# returns a generating set for the ideal of the ith exterior power of V_n
def ideal(n, i):
	# compute a generating set for the ideal
	genSet = []
	for j in range(n):
		for k in range(n):
			for l in range(n):
				if True: #(i<j and j<k) or (i<k and k<j) or (j<i and i<k):
					r = R(j, k, l)
					# we make use of the fact that x * R(j, k, l) is invariant 
					# under substituting (in x) elements appearing in R(j, k, l)
					#VnReduced = [v for v in V(n) if v != makePair(j, k) and v != makePair(k, l)]
					VnReduced = V(n)
					basisReduced = [ea.Element([(1, basisVector)]) for basisVector in ss.subsets(VnReduced, i - 2)]
					genSet.extend([x.mult(r) for x in basisReduced])

	# remove empty elements
	genSet = list(filter(lambda x: x.coeffVectorPairs, genSet))

	ea.getBasis(genSet)

	return genSet

# returns a matrix with columns that form a basis for the ideal of the ith exterior power of V_n
# passed the basis for the ideal and the dict of indices of subsets of the basis of V_n to save time
# result is in sparse CSR format
def idealBasisMatrix(n, i, idealBasis, indices):

	basisMatrix = lil_matrix((len(indices), len(idealBasis)))

	for j in range(len(idealBasis)):
		for k in range(len(idealBasis[j].coeffVectorPairs)):
			coeff, basisVector = idealBasis[j].coeffVectorPairs[k]
			basisMatrix[indices[tuple(basisVector)], j] = coeff

	return basisMatrix.tocsr()

# returns the image of idealBasisMatrix(n, i) under the action of perm
# passed the basis for the ideal and the dict of indices of subsets of the basis of V_n to save time
# result is in sparse CSC format
def permActionMatrix(n, i, idealBasis, indices, perm):
	actionMatrix = lil_matrix((len(indices), len(idealBasis)))

	for j in range(len(idealBasis)):
		imageVector = applyPermExtV(perm, idealBasis[j])
		for k in range(len(imageVector.coeffVectorPairs)):
			coeff, basisVector = imageVector.coeffVectorPairs[k]
			actionMatrix[indices[tuple(basisVector)], j] = coeff

	return actionMatrix.tocsc()

# returns a dict containing the values of the permutation character on V_n
def characterV(n):
	cycleTypes = Perm.partitions(1, n)
	charVVal = lambda x: (x.count(1) * (x.count(1) - 1)) / 2 + x.count(2)

	return {tuple(x): charVVal(x) for x in cycleTypes}

# returns a dict containing the values of the kth exterior power of the permutation character on V_n
def characterExtV(n, k):
	cycleTypes = Perm.partitions(1, n)

	if k == 0:
		return {tuple(cycleType): 1 for cycleType in cycleTypes}

	if k == 1:
		return characterV(n)

	previousChars = [characterExtV(n, i) for i in range(0,k)]
	characterExtVVal = lambda c: sum(map(lambda m: -sgn(m) * previousChars[1][tuple(Perm.powerCycleType(c, m))] * previousChars[k - m][tuple(c)], range(1,k + 1))) / k
	return {tuple(cycleType): characterExtVVal(cycleType) for cycleType in cycleTypes}

# returns the value of the the character of the ideal on a permutation
# passed the basis for the ideal and the dict of indices of subsets of the basis of V_n to save time
# also passed a factored version of idealBasisMatrix
def charVal(n, i, idealBasis, indices, iBMPseudoInv, perm):
	# if i < 2 the ideal is trivial
	if i < 2:
		return 0
	
	pAM = permActionMatrix(n, i, idealBasis, indices, perm)

	# compute the trace of inv(iBM)*permActionMatrix, which is the desired character
	# first handle case where iBMPseudoInv is 1xk
	if iBMPseudoInv.shape[0] == 1:
		trace = iBMPseudoInv.dot(pAM)
		if len(trace.shape) == 1:
			trace = trace[0]
		else:
			trace = trace[0,0]
	else:
		trace = 0
		for j in range(pAM.shape[1]):
			trace += iBMPseudoInv.getrow(j).dot(pAM.getcol(j))[0,0]

	return trace

# returns a list of character values for desired character on the cycle types in lexagraphical order
# the list output is so that it plays nicely with cp.interpolateCyclePolys
def character(n, i):
	cycleTypes = Perm.partitions(1, n)
	characterExtVDict = characterExtV(n, i)

	# if i < 2 the ideal is trivial
	if i < 2:
		return [characterExtVDict[tuple(cycleType)] for cycleType in cycleTypes]

	indices = ss.indices(V(n), i)

	print('calculating ideal basis')
	idealBasis = ideal(n, i)

	print('calculating pseudoinverse')
	iBM = idealBasisMatrix(n, i, idealBasis, indices)
	iBMNormal = iBM.T.dot(iBM)
	iBMPseudoInv = spsolve(iBMNormal, iBM.T)
	print('shape: {}, nonzero entries: {}'.format(iBMPseudoInv.shape, iBMPseudoInv.nnz))

	# define the character function
	char = lambda cycleType: characterExtVDict[tuple(cycleType)] - charVal(n, i, idealBasis, indices, iBMPseudoInv, Perm.fromCycleType(n, cycleType))

	print('calculating character')
	return parmap(char, cycleTypes)

# returns a dictionary from cycle types (as tuples) to character values for the desired character
def characterDict(n, i):
	return dict(zip(list(map(tuple, Perm.partitions(1, n)), character(n, i))))

# stores lists of character values in a given file for use later
# specifically, the characters on V_m through V_n
def characterDump(m, n, i, filename):
	output = open(filename, 'ab')
	for k in range(m, n + 1):
		characterValues = (k, i), character(k, i)
		print(characterValues)
		pickle.dump(characterValues, output)

	output.close()

# main function, runs characterDump
if __name__ == '__main__':
	characterDump(int(argv[1]), int(argv[2]), int(argv[3]), argv[4])
