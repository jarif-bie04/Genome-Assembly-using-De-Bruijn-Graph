from collections import defaultdict

class GenomeAssembler:
    def __init__(self, reads, k):
        self.reads = reads
        self.k = k

        self.kmers = []
        self.graph = defaultdict(set)

        self.in_degree = defaultdict(int)
        self.out_degree = defaultdict(int)

    # k-mers extraction
    def extract_kmers(self):
        self.kmers.clear()
        for read in self.reads:
            if len(read) < self.k:
                continue

            for i in range(len(read) - self.k + 1):
                self.kmers.append(read[i:i+self.k])

        return self.kmers

    # De Bruin Graph
    # De Bruijn Graph
    def build_graph(self):
        self.graph.clear()
        self.in_degree.clear()
        self.out_degree.clear()

        for kmer in self.kmers:
            left = kmer[:-1]
            right = kmer[1:]

            # Store only unique edges
            self.graph[left].add(right)

            self.out_degree[left] += 1
            self.in_degree[right] += 1

            if right not in self.graph:
                self.graph[right] = set()

        return self.graph

    # Start node finding
    def find_start_node(self):
        nodes = set(self.in_degree.keys()) | set(self.out_degree.keys())
        start = None
        for node in nodes:
            out_d = self.out_degree[node]
            in_d = self.in_degree[node]

            if out_d - in_d == 1:
                return node

            if out_d > 0:
                start = node

        return start

    # Eulerian path
    def eulerian_path(self):
        graph = {
            node: list(edges)
            for node, edges in self.graph.items()
        }
        start = self.find_start_node()
        if start is None:
            return []
        stack = [start]
        path = []
        while stack:
            node = stack[-1]
            if graph[node]:
                nxt = graph[node].pop()
                stack.append(nxt)
            else:
                path.append(stack.pop())
        path.reverse()

        return path

    # Genome Assemble
    def assemble_genome(self, path):
        if not path:
            return ""
        genome = path[0]
        for node in path[1:]:
            genome += node[-1]
        return genome

    # Statistics
    def statistics(self, path, genome):
        edges = sum(len(neighbors) for neighbors in self.graph.values())

        return {
            "reads": len(self.reads),
            "k": self.k,
            "nodes": len(self.graph),
            "edges": edges,
            "path_found": bool(path),
            "genome_length": len(genome)
        }

    def assemble(self):
        self.extract_kmers()
        self.build_graph()
        path = self.eulerian_path()
        genome = self.assemble_genome(path)
        stats = self.statistics(path, genome)
        return {
            "kmers": self.kmers,
            "graph": self.graph,
            "nodes": list(self.graph.keys()),
            "edges": [
                (u, v)
                for u in self.graph
                for v in self.graph[u]
            ],
            "path": path,
            "genome": genome,
            "stats": stats
        }