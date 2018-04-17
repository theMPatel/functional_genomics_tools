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