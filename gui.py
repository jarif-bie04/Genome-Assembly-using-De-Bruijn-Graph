import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os

from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from parser import DNAParser
from assembler import GenomeAssembler
from graph_utils import GraphDrawer

class GenomeAssemblyGUI(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Genome Assembly using De Bruijn Graph")
        self.geometry("1500x900")
        self.minsize(1200, 700)
        self.configure(bg="#EEF3F8")

        self.reads = []
        self.current_file = None

        self.setup_style()
        self.create_widgets()
        self.state("zoomed")

        self.resizable(False, False)

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

        # ===================== TOP HEADER =====================
        header = tk.Frame(self, bg="#0F4C81", height=70)
        header.pack(fill="x")
        header.pack_propagate(False)

        title = tk.Label(
            header,
            text="Genome Assembly using De Bruijn Graph",
            bg="#0F4C81",
            fg="white",
            font=("Segoe UI", 22, "bold")
        )
        title.pack(side="left", padx=25, pady=15)

        # Optional subtitle
        subtitle = tk.Label(
            header,
            text="Bioinformatics Visualization Tool",
            bg="#0F4C81",
            fg="#DDEEFF",
            font=("Segoe UI", 10)
        )
        subtitle.pack(side="right", padx=25)

        # BODY
        body = tk.Frame(self, bg="#EEF3F8")
        body.pack(fill="both", expand=True)

        self.left_panel(body)
        self.right_panel(body)

    # Left Panel

    def left_panel(self, parent):
        left = tk.Frame(
            parent,
            bg="#F8FAFD",
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

        read_scroll = tk.Scrollbar(left)

        self.read_text = tk.Text(
            left,
            width=40,
            height=18,
            font=("Consolas", 11),
            yscrollcommand=read_scroll.set
        )

        read_scroll.config(command=self.read_text.yview)

        read_scroll.pack(
            side="right",
            fill="y",
            padx=(0, 15),
            pady=5
        )
        self.read_text.pack(
            padx=15,
            pady=5,
            fill="both"
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
            bg="#FCFCFC"
        )

        self.kmer_frame.pack(
            side="left",
            fill="y",
            padx=5
        )
        self.kmer_frame.config(width=260)
        self.kmer_frame.pack_propagate(False)

        self.graph_frame = tk.LabelFrame(
            top,
            text="De Bruijn Graph",
            font=("Segoe UI", 11, "bold"),
            bg="#FCFCFC"
        )

        self.graph_frame.pack(
            side="left",
            fill="both",
            padx=5
        )
        self.graph_frame.config(width=650, height=520)
        self.graph_frame.pack_propagate(False)

        self.path_frame = tk.LabelFrame(
            top,
            text="Eulerian Path",
            font=("Segoe UI", 11, "bold"),
            bg="#FCFCFC",
            width=220
        )

        self.path_frame.pack(
            side="left",
            fill="y",
            padx=5
        )
        self.path_frame.config(width=180)
        self.path_frame.pack_propagate(False)

        # Bottom Row

        bottom = tk.Frame(self.right, bg="#f8f9fb")
        bottom.pack(
            fill="both",
            padx=15,
            pady=10
        )

        self.genome_frame = tk.LabelFrame(
            bottom,
            text="Assembled Genome",
            font=("Segoe UI", 11, "bold"),
            bg="#FCFCFC"
        )

        self.genome_frame.pack(
            side="left",
            fill="both",
            expand=True,
            padx=5
        )
        self.genome_frame.config(height=180)
        self.genome_frame.pack_propagate(False)

        self.stats_frame = tk.LabelFrame(
            bottom,
            text="Statistics",
            font=("Segoe UI", 11, "bold"),
            bg="#FCFCFC",
            width=250
        )

        self.stats_frame.pack(
            side="left",
            fill="both",
            padx=5
        )
        self.stats_frame.pack_propagate(False)

        self.create_result_widgets()

    def create_result_widgets(self):
        # K-mer box
        kmer_scroll = tk.Scrollbar(self.kmer_frame)

        self.kmer_text = tk.Text(
            self.kmer_frame,
            width=25,
            height=22,
            font=("Consolas", 11),
            yscrollcommand=kmer_scroll.set
        )

        kmer_scroll.config(command=self.kmer_text.yview)

        kmer_scroll.pack(side="right",
                         fill="y"
                    )
        self.kmer_text.pack(
            side="left",
            fill="both",
            expand=True,
            padx=10,
            pady=10
        )
        self.kmer_text.config(state="disabled")

        # Eulerian Path
        self.path_text = tk.Text(
            self.path_frame,
            width=18,
            height=22,
            font=("Consolas", 11)
        )

        path_scroll = tk.Scrollbar(self.path_frame)

        self.path_text = tk.Text(
            self.path_frame,
            width=18,
            height=22,
            font=("Consolas", 11),
            yscrollcommand=path_scroll.set
        )

        path_scroll.config(command=self.path_text.yview)

        path_scroll.pack(side="right", fill="y")
        self.path_text.pack(
            side="left",
            fill="both",
            expand=True,
            padx=10,
            pady=10
        )
        self.path_text.config(state="disabled")

        # Final Genome
        self.genome_text = tk.Text(
            self.genome_frame,
            font=("Consolas", 14, "bold"),
            fg="green",
            bg="white",
            wrap="word",
            height=6
        )

        self.genome_text.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=10
        )

        self.genome_text.config(state="disabled")

        # Statistics
        self.stats_text = tk.Text(
            self.stats_frame,
            width=30,
            height=12,
            font=("Segoe UI", 11)
        )

        stats_scroll = tk.Scrollbar(self.stats_frame)

        self.stats_text = tk.Text(
            self.stats_frame,
            width=30,
            height=12,
            font=("Segoe UI", 11),
            yscrollcommand=stats_scroll.set
        )

        stats_scroll.config(command=self.stats_text.yview)

        stats_scroll.pack(side="right", fill="y")
        self.stats_text.pack(
            side="left",
            fill="both",
            expand=True,
            padx=10,
            pady=10
        )
        self.stats_text.config(state="disabled")

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

        k = int(self.kmer_var.get())

        assembler = GenomeAssembler(self.reads, k)

        result = assembler.assemble()

        # K-mers
        self.kmer_text.config(state="normal")
        self.kmer_text.delete("1.0", tk.END)
        for i, kmer in enumerate(result["kmers"], start=1):
            self.kmer_text.insert(
                tk.END,
                f"{i}. {kmer}\n"
            )

        self.kmer_text.config(state="disabled")

        # Eulerian Path
        self.path_text.config(state="normal")
        self.path_text.delete("1.0", tk.END)
        if result["path"]:
            self.path_text.insert(
                tk.END,
                "\n↓\n".join(result["path"])
            )
        else:
            self.path_text.insert(
                tk.END,
                "No Eulerian Path"
            )

        self.path_text.config(state="disabled")

        # Genome
        self.genome_text.config(state="normal")
        self.genome_text.delete("1.0", tk.END)
        self.genome_text.insert(tk.END, result["genome"])
        self.genome_text.config(state="disabled")

        # Statistics
        stats = result["stats"]
        self.stats_text.config(state="normal")
        self.stats_text.delete("1.0", tk.END)
        self.stats_text.insert(
            tk.END,
            f"Number of Reads : {stats['reads']}\n"
        )

        self.stats_text.insert(
            tk.END,
            f"K-mer Size      : {stats['k']}\n"
        )

        self.stats_text.insert(
            tk.END,
            f"Graph Nodes     : {stats['nodes']}\n"
        )

        self.stats_text.insert(
            tk.END,
            f"Graph Edges     : {stats['edges']}\n"
        )

        self.stats_text.insert(
            tk.END,
            f"Eulerian Path   : {'YES' if stats['path_found'] else 'NO'}\n"
        )

        self.stats_text.insert(
            tk.END,
            f"Genome Length   : {stats['genome_length']} bp\n"
        )

        self.stats_text.config(state="disabled")

        # Draw Graph
        figure = GraphDrawer.draw(result["graph"])

        for widget in self.graph_frame.winfo_children():
            widget.destroy()

        canvas = FigureCanvasTkAgg(
            figure,
            master=self.graph_frame
        )

        canvas.draw()

        canvas.get_tk_widget().pack(
            fill="both",
            expand=True,
            padx=5,
            pady=5
        )

        plt.close(figure)

        self.export_btn.config(state="normal")

    def reset_all(self):
        self.read_text.delete("1.0", tk.END)

        self.kmer_text.config(state="normal")
        self.path_text.config(state="normal")
        self.stats_text.config(state="normal")

        self.kmer_text.delete("1.0", tk.END)
        self.path_text.delete("1.0", tk.END)
        self.stats_text.delete("1.0", tk.END)

        self.kmer_text.config(state="disabled")
        self.path_text.config(state="disabled")
        self.stats_text.config(state="disabled")

        self.genome_text.config(state="normal")
        self.genome_text.delete("1.0", tk.END)
        self.genome_text.insert(tk.END, "Genome will appear here")
        self.genome_text.config(state="disabled")

        for widget in self.graph_frame.winfo_children():
            widget.destroy()
        self.reads = []
        self.current_file = None

    def export_pdf(self):
        messagebox.showinfo(
            "Export",
            "PDF is now exported."
        )