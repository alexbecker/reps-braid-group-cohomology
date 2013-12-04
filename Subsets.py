# Subsets module
# codes functions for obtaining and indexing subsets of range(m, n)

# returns a list of all subsets of s of size k, in the form of sorted lists
# assumes s is sorted
def subsets(s, k):
	if len(s) < k:
		return []
	if k == 0:
		return [[]]

	smallSubsets = map(lambda x: [s[0]] + x, subsets(s[1:], k - 1))	# subsets starting with s[0]
	largeSubsets = subsets(s[1:], k)								# subsets starting with larger elts

	return smallSubsets + largeSubsets

# returns a dictionary from subsets of s of size k (as tuples) to indices
def indices(s, k):
	tuples = map(tuple, subsets(s, k))

	return dict(zip(tuples, range(len(tuples))))
