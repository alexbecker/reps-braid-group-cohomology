# Exterior Algebra module
# All methods are essentially sparse

# helper functions

# puts a basis vector in standard form, modifying the sign of the coefficient appropriately
# if the basis vector contains repeated terms (i.e. is 0), a coefficient of 0 is returns
# works by modified merge sort
def standardForm(coeffVectorPair):
	coeff, basisVector = coeffVectorPair
	length = len(basisVector)
	
	if length == 1:
		return coeffVectorPair
	
	firstHalf = basisVector[:(length+1)//2]
	secondHalf = basisVector[(length+1)//2:]
	sign1, firstHalfSorted = standardForm((1, firstHalf))
	sign2, secondHalfSorted = standardForm((1, secondHalf))

	# catch failure from firstHalf or secondHalf having repeated terms
	if sign1 == 0 or sign2 == 0:
		return 0, []

	# merge step
	coeff *= sign1 * sign2
	i = 0
	basisVectorSorted = []
	for elem in firstHalfSorted:
		# add elements from the second half until they are >= elem or we run out of elements
		while i < len(secondHalfSorted) and elem > secondHalfSorted[i]:
			basisVectorSorted.append(secondHalfSorted[i])
			i += 1
			coeff = -coeff

		# test for failure due to repeated elements
		if i < len(secondHalfSorted) and elem == secondHalfSorted[i]:
			return 0, []

		basisVectorSorted.append(elem)
	
	# append remaining elements to basisVectorSorted
	for elem in secondHalfSorted[i:]:
		basisVectorSorted.append(elem)
	
	return coeff, basisVectorSorted

# multiplies two basis vectors with coefficients
# assumes the arguments are in standard form
# result is also in standard form
# if result is 0, returns (0, [])
def mult(x, y):
	coeff1, basisVector1 = x
	coeff2, basisVector2 = y
	newBasisVector = []
	i = 0
	for elt in basisVector1:
		while i < len(basisVector2) and elt > basisVector2[i]:
			newBasisVector.append(basisVector2[i])
			i += 1
		if i < len(basisVector2) and elt == basisVector2[i]:
			return 0, []
		newBasisVector.append(elt)
	newBasisVector.extend(basisVector2[i:])

	if i % 2 == 0:
		return coeff1 * coeff2, newBasisVector
	else:
		return -coeff1 * coeff2, newBasisVector

# class Element
# represents an element of the exterior algebra of an ordered vector space

class Element:
	def __init__(self, coeffVectorPairs):
		self.coeffVectorPairs = coeffVectorPairs

	def __str__(self):
		return str(self.coeffVectorPairs)

	# standard form has all basis vectors in standard form, no repeated or null basis vectors, and the basis vectors sorted
	def standardForm(self, standardBasisVectors=False):
		# put the basis vectors in standard form, if not already
		if not standardBasisVectors:
			self.coeffVectorPairs = list(map(standardForm, self.coeffVectorPairs))
		
		# sort the result by the basis vector
		self.coeffVectorPairs.sort(key=lambda x: x[1])

		# add up elements with the same basis vector
		i = 0
		while i < len(self.coeffVectorPairs) - 1:
			coeff1, basisVector1 = self.coeffVectorPairs[i]
			coeff2, basisVector2 = self.coeffVectorPairs[i+1]
			if basisVector1 == basisVector2:
				self.coeffVectorPairs[i] = coeff1 + coeff2, basisVector1
				self.coeffVectorPairs.pop(i+1)
			else:
				i += 1

		# remove elements with coefficient 0
		self.coeffVectorPairs = [x for x in self.coeffVectorPairs if x[0] != 0]

	def scale(self, scalar):
		result = []
		for elem in self.coeffVectorPairs:
			coeff, basisVector = elem
			result.append((scalar * coeff, basisVector))

		return Element(result);

	# assumes both elements are in standard form
	# result is in standard form
	# if the result is 0, returns []
	def add(self, y):
		result = []
		i = 0
		for elem in self.coeffVectorPairs:
			coeff, basisVector = elem
			while i < len(y.coeffVectorPairs) and basisVector > y.coeffVectorPairs[i][1]:
				result.append(y.coeffVectorPairs[i])
				i += 1
			if i < len(y.coeffVectorPairs) and basisVector == y.coeffVectorPairs[i][1]:
				newCoeff = coeff + y.coeffVectorPairs[i][0]
				if newCoeff != 0:
					result.append((newCoeff, basisVector))
				i += 1
			else:
				result.append(elem)

		# append remaining element of y.coeffVectorPairs
		for elem in y.coeffVectorPairs[i:]:
			result.append(elem)

		return Element(result)

	# assumes both elements are in standard form (or at least their basisVectors are)
	# result is in standard form
	def mult(self, y):
		result = Element([mult(a, b) for a in self.coeffVectorPairs for b in y.coeffVectorPairs])
		result.standardForm(standardBasisVectors = True)

		return result

	# returns the coefficient of a given basis vector in an element
	def getCoeff(self, basisVector):
		for elem in self.coeffVectorPairs:
			coeff, vector = elem
			if vector == basisVector:
				return coeff

		return 0

# other useful functions

# reduces a list of elements to a basis for their span
def getBasis(spanningSet):
	i = 0
	while i < len(spanningSet):
		# update user on progress
		if i % 100 == 0:
			print('{} / {}'.format(i, len(spanningSet)))

		# clear the column corresponding to the first basis vector of spanningSet[i]
		coeff, basisVector = spanningSet[i].coeffVectorPairs[0]
		j = i + 1
		while j < len(spanningSet):
			multiplier = - spanningSet[j].getCoeff(basisVector) / float(coeff)
			if multiplier != 0:
				spanningSet[j] = spanningSet[j].add(spanningSet[i].scale(multiplier))
		
			#remove elements which are 0
			if not spanningSet[j].coeffVectorPairs:
				del spanningSet[j]
			else:
				j += 1
		i += 1
