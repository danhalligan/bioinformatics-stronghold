from collections import defaultdict


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
