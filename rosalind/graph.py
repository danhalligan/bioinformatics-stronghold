import re

from collections import defaultdict
from itertools import permutations


class Graph:
    """Store a graph in a dictionary where key is node and value is a list of
    connections"""

    def __init__(self, adjacency_list):
        self.graph = defaultdict(list)
        for x in adjacency_list:
            self.graph[x[0]].append(x[1])
            self.graph[x[1]].append(x[0])
        self.nodes = list(self.graph.keys())

    def count_distinct(self):
        visited = {}
        n_graphs = 0

        def visit_nodes(node):
            visited[node] = True
            for i in self.graph[node]:
                if i not in visited:
                    visit_nodes(i)

        for node in list(self.nodes):
            if node not in visited:
                visit_nodes(node)
                n_graphs += 1

        return n_graphs


def overlap_graph(seqs, n=3):
    """Build an overlap graph of sequences where prefix of A matches suffix of
    B"""
    for pair in permutations(seqs, 2):
        if pair[0].seq.endswith(pair[1].seq[:n]):
            yield (pair[0].id, pair[1].id)


class Node:
    """Store a node in a graph, retaining parent and child relationships"""

    def __init__(self, id, name, parent, children=[]):
        self.id = id
        self.name = name
        self.parent = parent
        self.children = children

    # Method to get the depth of the node (for printing)
    def depth(self):
        current_node = self
        depth = 0
        while type(current_node.parent):
            current_node = current_node.parent
            depth += 1
        return depth

    def __repr__(self):
        if len(self.children) > 0:
            children = ",".join(x.__repr__() for x in self.children)
            return "({}){}".format(children, self.name)
        else:
            return self.name


def parse(newick):
    # find a name or a delimiter.
    # a name is any A-z_ characters followed by [,;)]
    tokens = re.finditer(r"([A-z_]*)([,;)])|(\S)", newick)

    def recurse(nextid, pid):
        id = nextid
        children = []

        name, delim, ch = next(tokens).groups()
        while ch == "(":
            while ch in "(,":
                node, ch, nextid = recurse(nextid + 1, id)
                children.append(node)
            name, delim, ch = next(tokens).groups()
        return Node(id, name, pid, children), delim, nextid

    return recurse(0, -1)[0]


newick = "(long_name,Basd,(C,D)E,G),(H,I);"
parse(newick)

# def parse(newick):
#     def recurse(nextid, pid, i):
#         id = nextid;
#         children = []
#         ch = newick[i]
#         m = re.match(r"^[A-z_]+", newick[i+1:])
#         if m:
#             name = m.string
#         if ch == "(":
#             while ch in "(,":
#                 node, ch, nextid = recurse(nextid+1, id)
#                 children.append(node)
#             name, delim, ch = next(tokens).groups()
#         else:

#         return Node(id, name, pid, children), delim, nextid

#     return recurse(0, -1)[0]


# def nwck(tree, nodes):
#     """Distances in Trees"""
#     i = 0
#     node = 0
#     graph = {}
#     while i < len(tree):
#         if tree[i] == "(":
#             print(i)
#         i += 1


# class Node:
#     def __init__(self, name):
#         self.name = name
#         self.children = []
#         self.parent = None

#     # Method to get the depth of the node (for printing)
#     def get_depth(self):
#         current_node = self
#         depth = 0
#         while current_node.parent:
#             current_node = current_node.parent
#             depth += 1
#         return depth

#     # String representation
#     def __str__(self):
#         return self.name

# newick = "(A,Bcad,(C,D)E,G)F"

# root = None
# na = ""
# stack = []
# for i in list(reversed(newick)):
#     if i == ')':
#         if na != "":
#             node = Node(na)
#             na = ""
#             if len(stack):
#                 stack[-1].children.append(node)
#                 node.parent = stack[-1]
#             else:
#                 root = node
#             stack.append(node)
#     elif i == '(':
#         if (na != ""):
#             node = Node(na)
#             na = ""
#             stack[-1].children.append(node)
#             node.parent = stack[-1]
#         stack.pop()
#     elif i == ',':
#         if (na != ""):
#             node = Node(na)
#             na = ""
#             stack[-1].children.append(node)
#             node.parent = stack[-1]
#     else:
#         na += i

# # Just to print the parsed tree.
# print_stack = [root]
# while len(print_stack):
#     node = print_stack.pop()
#     print("-" * node.get_depth(), node)
#     print_stack.extend(node.children)
