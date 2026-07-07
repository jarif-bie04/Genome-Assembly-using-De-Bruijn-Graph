import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os

from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas

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

        # Labels
        style.configure(
            "Title.TLabel",
            font=("Segoe UI", 18, "bold"),
            foreground="#1F4E79",
            background="#EEF3F8"
        )

        style.configure(
            "Header.TLabel",
            font=("Segoe UI", 12, "bold"),
            foreground="#1F4E79",
            background="#EEF3F8"
        )

        style.configure(
            "Normal.TLabel",
            font=("Segoe UI", 10),
            background="#EEF3F8"
        )

        # Blue Button
        style.configure(
            "Blue.TButton",
            font=("Segoe UI", 10, "bold"),
            foreground="white",
            background="#1976D2",
            borderwidth=0,
            padding=10
        )

        style.map(
            "Blue.TButton",
            background=[
                ("active", "#1565C0"),
                ("pressed", "#0D47A1")
            ],
            foreground=[
                ("disabled", "#CCCCCC")
            ]
        )

        # Green Button
        style.configure(
            "Green.TButton",
            font=("Segoe UI", 10, "bold"),
            foreground="white",
            background="#2E7D32",
            borderwidth=0,
            padding=10
        )

        style.map(
            "Green.TButton",
            background=[
                ("active", "#1B5E20"),
                ("pressed", "#0B3D16")
            ]
        )

        # Red Button
        style.configure(
            "Red.TButton",
            font=("Segoe UI", 10, "bold"),
            foreground="white",
            background="#D32F2F",
            borderwidth=0,
            padding=10
        )

        style.map(
            "Red.TButton",
            background=[
                ("active", "#C62828"),
                ("pressed", "#8E0000")
            ]
        )

        # Combobox
        style.configure(
            "TCombobox",
            padding=6,
            font=("Segoe UI", 10)
        )

    # Main Panel

    def create_widgets(self):

        # TOP HEADER
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

        # BODY
        body = tk.Frame(self, bg="#EEF3F8")
        body.pack(fill="both", expand=True)

        self.left_panel(body)
        self.right_panel(body)

        # FOOTER
        footer = tk.Frame(self, bg="#0F4C81", height=30)
        footer.pack(fill="x", side="bottom")
        footer.pack_propagate(False)

        tk.Label(
            footer,
            text="Developed by: Jarif Ayman © 2026 Bioinformatics Engineering, BAU",
            bg="#0F4C81",
            fg="white",
            font=("Segoe UI", 9)
        ).pack(expand=True)

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

        top = tk.Frame(self.right, bg="#EEF3F8")
        top.pack(fill="both", padx=15, pady=10)

        top.grid_columnconfigure(0, weight=1)
        top.grid_columnconfigure(1, weight=3)
        top.grid_columnconfigure(2, weight=1)

        top.pack(
            fill="both",
            expand=True,
            padx=15
        )

        self.kmer_frame = tk.Frame(
            top,
            bg="white",
            bd=0,
            highlightbackground="#D9E3F0",
            highlightthickness=1
        )

        tk.Label(
            self.kmer_frame,
            text="K-mers Extracted",
            bg="white",
            fg="#214A9A",
            font=("Segoe UI", 11, "bold")
        ).pack(anchor="w", padx=12, pady=(8, 4))

        self.kmer_frame.grid(
            row=0,
            column=0,
            sticky="nsew",
            padx=5
        )

        self.kmer_frame.config(width=260)
        self.kmer_frame.pack_propagate(False)

        graph_card = tk.Frame(
            top,
            bg="white",
            bd=1,
            relief="solid"
        )

        graph_card.grid(
            row=0,
            column=1,
            sticky="nsew",
            padx=8
        )

        graph_card.grid_propagate(False)
        graph_card.config(width=700, height=520)

        tk.Label(
            graph_card,
            text="De Bruijn Graph ((k-1)-mers)",
            bg="white",
            fg="#1A3E8E",
            font=("Segoe UI", 12, "bold")
        ).pack(pady=10)

        self.graph_frame = tk.Frame(
            graph_card,
            bg="white"
        )

        self.graph_frame.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=(0, 10)
        )

        path_card = tk.Frame(
            top,
            bg="white",
            bd=1,
            relief="solid"
        )

        path_card.grid(
            row=0,
            column=2,
            sticky="nsew",
            padx=8
        )

        path_card.grid_propagate(False)
        path_card.config(width=220, height=520)

        tk.Label(
            path_card,
            text="Eulerian Path",
            bg="white",
            fg="#1A3E8E",
            font=("Segoe UI", 12, "bold")
        ).pack(pady=10)

        self.path_frame = tk.Frame(
            path_card,
            bg="white"
        )

        self.path_frame.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=(0, 10)
        )

        # Bottom Row

        bottom = tk.Frame(self.right, bg="#f8f9fb")
        bottom.pack(
            fill="both",
            padx=15,
            pady=10
        )

        genome_card = tk.Frame(
            bottom,
            bg="white",
            bd=1,
            relief="solid"
        )

        genome_card.pack(
            side="left",
            fill="both",
            expand=True,
            padx=8
        )

        genome_card.pack_propagate(False)
        genome_card.config(height=200)

        tk.Label(
            genome_card,
            text="Assembled Genome (Final Result)",
            bg="white",
            fg="#1A3E8E",
            font=("Segoe UI", 12, "bold")
        ).pack(anchor="w", padx=15, pady=10)

        self.genome_frame = tk.Frame(
            genome_card,
            bg="white"
        )

        self.genome_frame.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=(0, 10)
        )

        stats_card = tk.Frame(
            bottom,
            bg="white",
            bd=1,
            relief="solid"
        )

        stats_card.pack(
            side="left",
            fill="both",
            padx=8
        )

        stats_card.pack_propagate(False)
        stats_card.config(width=300, height=200)

        tk.Label(
            stats_card,
            text="Assembly Statistics",
            bg="white",
            fg="#1A3E8E",
            font=("Segoe UI", 12, "bold")
        ).pack(anchor="w", padx=15, pady=10)

        self.stats_frame = tk.Frame(
            stats_card,
            bg="white"
        )

        self.stats_frame.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=(0, 10)
        )

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
        print("Export Enabled")

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

        filename = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF Files", "*.pdf")],
            title="Save PDF"
        )

        if not filename:
            return

        c = canvas.Canvas(filename, pagesize=A4)

        y = 800

        c.setFont("Helvetica-Bold", 16)
        c.drawString(50, y, "Genome Assembly Report")
        y -= 40

        c.setFont("Helvetica", 11)

        c.drawString(50, y, "Input Reads:")
        y -= 20

        for read in self.reads:
            c.drawString(70, y, read)
            y -= 15

            if y < 60:
                c.showPage()
                y = 800

        y -= 20

        c.setFont("Helvetica-Bold", 12)
        c.drawString(50, y, "Assembled Genome")
        y -= 20

        genome = self.genome_text.get("1.0", tk.END).strip()

        c.setFont("Courier", 10)
        c.drawString(70, y, genome)
        y -= 30

        c.setFont("Helvetica-Bold", 12)
        c.drawString(50, y, "Statistics")
        y -= 20

        stats = self.stats_text.get("1.0", tk.END).strip()

        c.setFont("Courier", 10)

        for line in stats.split("\n"):
            c.drawString(70, y, line)
            y -= 15

        c.save()

        messagebox.showinfo(
            "Success",
            f"PDF saved successfully.\n\n{filename}"
        )