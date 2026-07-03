from assembler import GenomeAssembler
from graph_utils import GraphDrawer

reads = [
    "ATGCG",
    "TGCGA",
    "GCGAT",
    "CGATC",
    "GATCG",
    "ATCGA"
]

k = 4

assembler = GenomeAssembler(reads, k)

result = assembler.assemble()

print("K-mers:")
print(result["kmers"])

print("\nGenome:")
print(result["genome"])

print("\nStatistics:")
print(result["stats"])

GraphDrawer.draw(result["graph"])