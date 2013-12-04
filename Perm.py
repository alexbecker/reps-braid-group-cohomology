# Perm module
# implements functions relating to permutations, which are stored as lists

from fractions import gcd

# given a list of cycle lengths, returns a permutation of range(n) of the given cycle type
def fromCycleType(n, cycleType):
	result = []
	i = 0
	for length in cycleType:
		result.extend(range(i + 1, i + length) + [i])
		i += length

	return result

# returns a list of all partitions of n with smallest element at least m
def partitions(m, n):
	if n == 0:
		return [[]]

	return [[k] + partition for k in range(m, n + 1) for partition in partitions(k, n - k)]

# returns a list of permutations of range(n) containing one of each cycle type
def allCycleTypes(n):
	return map(lambda x: fromCycleType(n, x), partitions(1, n))

# returns the cycle type of the kth power of permutation of a given cycle type
def powerCycleType(cycleType, k):
	result = []
	for length in cycleType:
		result.extend([length / gcd(length, k) for i in range(gcd(length, k))])

	result.sort()

	return result
