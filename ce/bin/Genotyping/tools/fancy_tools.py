###################################################################
#
# Fancy tools!
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

class Disjointset(object):
    # Disjoint sets object. 
    # Balances by:
    #   -> Rank heuristic
    #   -> Path compression

    def __init__(self, size):
        self._size = size
        self._rank = [0] * self._size
        self._parents = list(range(size))

    def get_parent(self, index):
        # This will flatten the disjoint set so that all
        # children point to the true parent
        if self._parents[index] != index:
            self._parents[index] = self.get_parent(self._parents[index])

        return self._parents[index]

    def merge(self, source, destination):
        # Merges based on rank heuristic
        real_source = self.get_parent(source)
        real_destination = self.get_parent(destination)

        if real_source == real_destination:
            return

        if self._rank[real_destination] > self._rank[real_source]:
            self._parents[real_source] = real_destination

        else:
            self._parents[real_destination] = real_source

            if self._rank[real_destination] == self._rank[real_source]:
                self._rank[real_source] += 1

class DecisionTree(object):

    def __init__(self):
        self._attr = None
        self._value = None
        self._tree = {}

    def update(self, branch):
        # Attrs should come in transposed and sorted by column header
        # values of attrs can only be True or False

        if not self._attr:
            self._attr = branch[0][0]

        else:
            assert self._attr == branch[0][0]

        if isinstance(branch[0][1], basestring):
            self._value = branch[0][1]

        elif branch[0][1] in self._tree:
            self._tree[branch[0][1]].update(branch[1:])

        else:
            self._tree[branch[0][1]] = DecisionTree()
            self._tree[branch[0][1]].update(branch[1:])

    def recurse(self, branch):

        if self._value:
            return self._value

        elif not branch[0][0] == self._attr:
            return

        else:
            return self._tree[branch[0][1]].recurse(branch[1:])

def edit_distance(c, n, matrix=False):
    c_l = len(c)
    n_l = len(n)

    if not c_l:
        return n_l

    if not n_l:
        return c_l
        
    if c_l < n_l:
        return edit_distance(n, c)

    m = [[0]*(n_l+1) for _ in range(c_l+1)]

    for i in range(c_l+1):
        m[i][0] = i

    for j in range(n_l+1):
        m[0][j] = j

    for i in range(c_l):

        for j in range(n_l):

            if c[i] == n[j]:

                m[i+1][j+1] = m[i][j]

            else:

                m[i+1][j+1] = min(
                    m[i][j],
                    m[i+1][j],
                    m[i][j+1]
                ) + 1

    if matrix:
        return m
    else:
        return m[c_l][n_l]

def build_consensus(c, n):
    # Get the edit distance_matrix 
    if not c:
        return n

    if not n:
        return c

    m = edit_distance(c, n, matrix=True)
    
    # Here we will rebuild the sequence
    i = len(m)
    j = len(m[0])

    final = []
    while i >= 1 and j >= 1:

        v = m[i][j]

        up_left = m[i-1][j-1]
        left = m[i-1][j]
        up = m[i][j-1]

        if all(up_left == x for x in [left, up]):
            i-=1
            j-=1
            final.append(n[j])
        
        else:
            m_val = min([up_left, left, up])

            if up_left == m_val:

                i -= 1
                j -= 1

                if c[i] == n[j]:
                    final.append(n[j])
                else:
                    final.append('N')

            elif left == m_val:

                i-=1
                final.append('N')

            else:
                j -= 1
                final.append('N')

    while i >= 1:
        final.append('N')
        i-=1

    while j >= 1:
        final.append('N')
        j-=1

    final.reverse()

    return final
