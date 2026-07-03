import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os

from parser import DNAParser


class GenomeAssemblyGUI(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Genome Assembly using De Bruijn Graph")
        self.geometry("1500x900")
        self.minsize(1200, 700)
        self.configure(bg = "#f2f4f7")

        self.reads = []
        self.current_file = None

        self.setup_style()
        self.create_widgets()

    def setup_style(self):
        style = ttk.Style()
        style.theme_use("clam")

        style.configure(
            "Title.TLabel",
            font=("Segoe UI", 18, "bold"),
            foreground="#1f4e79",
            background="#f2f4f7"
        )

        style.configure(
            "Header.TLabel",
            font=("Segoe UI", 12, "bold"),
            foreground="#1f4e79",
            background="#f2f4f7"
        )

        style.configure(
            "Normal.TLabel",
            font=("Segoe UI", 10),
            background="#f2f4f7"
        )

        style.configure(
            "Blue.TButton",
            font=("Segoe UI", 10, "bold"),
            padding=8
        )

        style.configure(
            "Green.TButton",
            font=("Segoe UI", 10, "bold"),
            padding=8
        )

        style.configure(
            "Red.TButton",
            font=("Segoe UI", 10, "bold"),
            padding=8
        )

    # Main Panel

    def create_widgets(self):
        title = ttk.Label(
            self,
            text="Genome Assembly using De Bruijn Graph",
            style="Title.TLabel"
        )

        title.pack(pady=10)
        body = tk.Frame(self, bg="#f2f4f7")
        body.pack(fill="both", expand=True)
        self.left_panel(body)
        self.right_panel(body)

    # Left Panel

    def left_panel(self, parent):
        left = tk.Frame(
            parent,
            bg="white",
            bd=1,
            relief="solid",
            width=380
        )

        left.pack(side="left", fill="y", padx=10, pady=10)

        # Input field

        ttk.Label(
            left,
            text="INPUT DNA READS",
            style="Header.TLabel"
        ).pack(anchor="w", padx=15, pady=(15, 5))

        self.read_text = tk.Text(
            left,
            width=40,
            height=18,
            font=("Consolas", 11)
        )

        self.read_text.pack(
            padx=15,
            pady=5,
            fill="x"
        )

        # load Button

        upload_btn = ttk.Button(
            left,
            text="Load FASTA / FASTQ",
            command=self.load_file,
            style="Blue.TButton"
        )

        upload_btn.pack(
            padx=15,
            pady=10,
            fill="x"
        )

        # k-mer dropdown

        ttk.Label(
            left,
            text="SELECT K-MER SIZE",
            style="Header.TLabel"
        ).pack(anchor="w", padx=15, pady=(20, 5))

        self.kmer_var = tk.StringVar()

        self.kmer_combo = ttk.Combobox(
            left,
            textvariable=self.kmer_var,
            state="readonly",
            values=[3, 4, 5, 6, 7, 8, 9, 10]
        )

        self.kmer_combo.current(1)

        self.kmer_combo.pack(
            padx=15,
            fill="x"
        )

        ttk.Label(
            left,
            text="ACTIONS",
            style="Header.TLabel"
        ).pack(anchor="w", padx=15, pady=(25, 5))

        # Show Result button
        self.show_btn = ttk.Button(
            left,
            text="Show Result",
            command=self.show_result,
            style="Blue.TButton"
        )

        self.show_btn.pack(
            padx=15,
            pady=5,
            fill="x"
        )

        # Reset Button
        self.reset_btn = ttk.Button(
            left,
            text="Reset",
            command=self.reset_all,
            style="Green.TButton"
        )

        self.reset_btn.pack(
            padx=15,
            pady=5,
            fill="x"
        )

        # Export Button
        self.export_btn = ttk.Button(
            left,
            text="Export PDF",
            state="disabled",
            command=self.export_pdf,
            style="Red.TButton"
        )

        self.export_btn.pack(
            padx=15,
            pady=5,
            fill="x"
        )

    # Right Panel
    def right_panel(self, parent):
        self.right = tk.Frame(
            parent,
            bg="#f8f9fb",
            bd=1,
            relief="solid"
        )

        self.right.pack(
            side="left",
            fill="both",
            expand=True,
            padx=10,
            pady=10
        )

        ttk.Label(
            self.right,
            text="RESULTS",
            style="Header.TLabel"
        ).pack(anchor="w", padx=20, pady=(15, 10))

        # TOP ROW

        top = tk.Frame(self.right, bg="#f8f9fb")
        top.pack(fill="both", expand=True, padx=15)

        self.kmer_frame = tk.LabelFrame(
            top,
            text="K-mers Extracted",
            font=("Segoe UI", 11, "bold"),
            bg="white"
        )

        self.kmer_frame.pack(side="left", fill="both", expand=True, padx=5)

        self.graph_frame = tk.LabelFrame(
            top,
            text="De Bruijn Graph",
            font=("Segoe UI", 11, "bold"),
            bg="white"
        )

        self.graph_frame.pack(side="left", fill="both", expand=True, padx=5)

        self.path_frame = tk.LabelFrame(
            top,
            text="Eulerian Path",
            font=("Segoe UI", 11, "bold"),
            bg="white",
            width=220
        )

        self.path_frame.pack(side="left", fill="y", padx=5)

        # Bottom Row

        bottom = tk.Frame(self.right, bg="#f8f9fb")
        bottom.pack(fill="x", padx=15, pady=15)

        self.genome_frame = tk.LabelFrame(
            bottom,
            text="Assembled Genome",
            font=("Segoe UI", 11, "bold"),
            bg="white"
        )

        self.genome_frame.pack(
            side="left",
            fill="both",
            expand=True,
            padx=5
        )

        self.stats_frame = tk.LabelFrame(
            bottom,
            text="Statistics",
            font=("Segoe UI", 11, "bold"),
            bg="white",
            width=300
        )

        self.stats_frame.pack(
            side="left",
            fill="y",
            padx=5
        )

        self.create_result_widgets()

    def create_result_widgets(self):
        # K-mer box
        self.kmer_text = tk.Text(
            self.kmer_frame,
            width=25,
            height=22,
            font=("Consolas", 11)
        )

        self.kmer_text.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=10
        )

        # Graph
        self.graph_placeholder = tk.Label(
            self.graph_frame,
            text="\n\nGraph will appear here",
            font=("Segoe UI", 13),
            bg="white",
            fg="gray"
        )

        self.graph_placeholder.pack(expand=True)

        # Eulerian Path
        self.path_text = tk.Text(
            self.path_frame,
            width=18,
            height=22,
            font=("Consolas", 11)
        )

        self.path_text.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=10
        )

        # Final Genome
        self.genome_label = tk.Label(
            self.genome_frame,
            text="Genome will appear here",
            font=("Consolas", 24, "bold"),
            bg="white",
            fg="green"
        )

        self.genome_label.pack(
            expand=True,
            pady=30
        )

        # Statistics
        self.stats_text = tk.Text(
            self.stats_frame,
            width=30,
            height=12,
            font=("Segoe UI", 11)
        )

        self.stats_text.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=10
        )

    def load_file(self):
        filepath = filedialog.askopenfilename(
            title="Open DNA File",
            filetypes=[
                ("FASTA Files", "*.fa *.fasta"),
                ("FASTQ Files", "*.fq *.fastq"),
                ("All Files", "*.*")
            ]
        )

        if not filepath:
            return

        try:
            self.reads = DNAParser.load(filepath)
            self.current_file = filepath
            self.read_text.delete("1.0", tk.END)
            for read in self.reads:
                self.read_text.insert(tk.END, read + "\n")
            messagebox.showinfo(
                "Loaded",
                f"{len(self.reads)} reads loaded successfully."
            )

        except Exception as e:
            messagebox.showerror(
                "Error",
                str(e)
            )

    def show_result(self):
        manual_reads = DNAParser.parse_manual(
            self.read_text.get("1.0", tk.END)
        )
        if manual_reads:
            self.reads = manual_reads
        if not self.reads:
            messagebox.showwarning(
                "Warning",
                "Please enter DNA reads."
            )
            return
        self.kmer_text.delete("1.0", tk.END)
        k = int(self.kmer_var.get())
        total = 0
        for read in self.reads:
            for i in range(len(read) - k + 1):
                kmer = read[i:i + k]
                total += 1
                self.kmer_text.insert(
                    tk.END,
                    f"{total}. {kmer}\n"
                )
        self.stats_text.delete("1.0", tk.END)
        self.stats_text.insert(
            tk.END,
            f"Number of Reads : {len(self.reads)}\n"
        )
        self.stats_text.insert(
            tk.END,
            f"K-mer Size      : {k}\n"
        )
        self.stats_text.insert(
            tk.END,
            "Graph Nodes     : --\n"
        )
        self.stats_text.insert(
            tk.END,
            "Graph Edges     : --\n"
        )
        self.stats_text.insert(
            tk.END,
            "Eulerian Path   : --\n"
        )
        self.stats_text.insert(
            tk.END,
            "Genome Length   : --\n"
        )

    def reset_all(self):
        self.read_text.delete("1.0", tk.END)
        self.kmer_text.delete("1.0", tk.END)
        self.path_text.delete("1.0", tk.END)
        self.stats_text.delete("1.0", tk.END)
        self.genome_label.config(
            text="Genome will appear here"
        )
        self.reads = []
        self.current_file = None

    def export_pdf(self):
        messagebox.showinfo(
            "Export",
            "PDF is now exported."
        )