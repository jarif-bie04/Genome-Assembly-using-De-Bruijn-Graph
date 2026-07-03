import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.pyplot import arrow


class GraphDrawer:

    @staticmethod
    def draw(graph, save_path=None):
        G = nx.DiGraph()

        for source in graph:
            G.add_node(source)

            for target in graph[source]:
                G.add_edge(source, target)

        plt.figure(figsize=(9, 7))

        # Layout
        pos = nx.spring_layout(
            G,
            seed=42,
            k=3.5,
            iterations=300
        )

        # Nodes
        nx.draw_networkx_nodes(
            G,
            pos,
            node_size=600,
            node_color="#cfe8ff",
            edgecolors="#1f77b4",
            linewidths=2
        )

        # Labels
        nx.draw_networkx_labels(
            G,
            pos,
            font_size=9,
            font_weight="bold"
        )

        # Directed edges
        ax = plt.gca()

        for u, v in G.edges():
            ax.annotate(
                "",
                xy=pos[v],
                xytext=pos[u],
                arrowprops=dict(
                    arrowstyle="-|>",
                    color="black",
                    lw=2,
                    mutation_scale=20,
                    shrinkA=20,
                    shrinkB=20,
                ),
            )

        plt.title(
            "De Bruijn Graph",
            fontsize=15,
            weight="bold"
        )

        plt.axis('off')

        if save_path:
            plt.savefig(
                save_path,
                dpi=300,
                bbox_inches="tight"
            )

        figure = plt.gcf()
        return figure