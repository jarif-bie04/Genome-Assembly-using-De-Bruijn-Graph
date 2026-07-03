import os

class DNAParser:

    VALID = {'A', 'T', 'C', 'G'}

    @staticmethod
    def clean_sequence(seq):
        seq = seq.upper()
        seq = seq.replace(" ", "")
        seq = seq.replace("\n", "")
        seq = seq.replace("\r", "")

        return seq

    @staticmethod
    def validate_sequence(seq):
        return all(base in DNAParser.VALID for base in seq)

    @staticmethod
    def parse_manual(text):
        reads = []
        lines = text.strip().splitlines()

        for line in lines:
            line = DNAParser.clean_sequence(line)
            if line == "":
                continue
            if DNAParser.validate_sequence(line):
                reads.append(line)
        return reads

    @staticmethod
    def parse_fasta(filepath):
        reads = []
        current = ""
        with open(filepath, "r") as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current:
                        current = DNAParser.clean_sequence(current)
                        if DNAParser.validate_sequence(current):
                            reads.append(current)
                    current = ""
                else:
                    current += line
        if current:
            current = DNAParser.clean_sequence(current)
            if DNAParser.validate_sequence(current):
                reads.append(current)
        return reads

    @staticmethod
    def parse_fastq(filepath):
        reads = []
        with open(filepath, "r") as file:
            lines = file.readlines()
        for i in range(1, len(lines), 4):
            seq = DNAParser.clean_sequence(lines[i])
            if DNAParser.validate_sequence(seq):
                reads.append(seq)
        return reads

    @staticmethod
    def load(filepath):
        ext = os.path.splitext(filepath)[1].lower()
        if ext in [".fa", ".fasta"]:
            return DNAParser.parse_fasta(filepath)
        elif ext in [".fq", ".fastq"]:
            return DNAParser.parse_fastq(filepath)
        else:
            raise ValueError("Unsupported file format.")