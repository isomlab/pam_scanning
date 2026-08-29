"""Tkinter GUI front-end for PAM-scanning.

Collects run parameters through a themed form with hover help and hands them to
:func:`pam_scanning.chimeras.pamscan`. One or more ORFs can be queued; each ORF
has its own gene name and ORF FASTA, plus either its own 5'/3' flanks (per-ORF
mode) or a single global 5'/3' flank pair shared by every ORF (global mode). Each
global flank is taken from a FASTA file or from a sequence typed/pasted straight
into the form. The genome, codon table, primer suffixes, and scan parameters are
always shared. ORFs can be added one at a time or discovered from a folder.
Launch with ``pam-scan-gui``.
"""

import contextlib
import json
import os
import queue
import threading
import time
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, font as tkfont
from pathlib import Path

from pam_scanning import __version__, fetch_cds
from pam_scanning.chimeras import parse_sequence_text, parse_codon_positions
from pam_scanning.cli import blast_db_prefix, discover_orf_folder, gene_name_from_orf_path

# FASTA extensions the folder loader inspects (matches cli.discover_orf_folder).
_FASTA_EXTS = (".fa", ".fasta", ".fna", ".fas")


class _ConsoleWriter:
    """A minimal stdout-like sink that forwards written text to a callback.

    Used with :func:`contextlib.redirect_stdout` so the pipeline's ``print``
    progress is captured and shown in the GUI's console instead of the terminal.
    The callback is expected to be thread-safe (it enqueues for the UI thread).
    """

    def __init__(self, sink):
        self._sink = sink

    def write(self, text):
        if text:
            self._sink(text)
        return len(text)

    def flush(self):
        pass

# Where the last-used dialog directory is remembered between sessions.
_STATE_PATH = Path.home() / ".pam_scanning" / "gui_state.json"


def _load_last_dir():
    """Return the directory the last file dialog used, or the cwd if none/stale."""
    try:
        saved = json.loads(_STATE_PATH.read_text()).get("last_dir")
    except (OSError, ValueError):
        saved = None
    return saved if saved and os.path.isdir(saved) else os.getcwd()


def _save_last_dir(path):
    """Persist the last-used dialog directory; ignore any I/O failure."""
    try:
        _STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
        _STATE_PATH.write_text(json.dumps({"last_dir": path}))
    except OSError:
        pass

# ---------------------------------------------------------------------------
# Palette and field definitions
# ---------------------------------------------------------------------------

BG = "#eef2f5"          # window background
CARD = "#ffffff"        # section card background
HEADER_BG = "#1f3a5f"   # header banner
HEADER_FG = "#ffffff"
SUBTITLE_FG = "#c3d0de"
TEXT = "#1f2a36"
MUTED = "#6b7a8d"
ACCENT = "#2e7d32"      # primary action (Run)
ACCENT_ACTIVE = "#256628"
TIP_BG = "#fffbe6"
TIP_BORDER = "#c9b458"
INVALID = "#b00020"     # sequence box holds non-DNA characters
# Scrolling (Tk 9 delivers trackpad scroll as <TouchpadScroll>, decoded to pixel deltas).
TOUCHPAD_SCALE = 1.4      # pixels scrolled per unit of trackpad delta
MOUSE_WHEEL_PIXELS = 40   # pixels per mouse-wheel notch
CONSOLE_BG = "#0f1b28"  # progress console background
CONSOLE_FG = "#d7e0ea"  # progress console text
CONSOLE_MUTED = "#7f93a8"
PLACEHOLDER = "No file selected"
# Shown by the codon picker while nothing is chosen; the "(optional)" matches the
# file rows so both codon-selection controls read the same way.
NO_CODONS_PICKED = "No codons picked graphically  (optional)"

# Per-ORF flank inputs (used only in per-ORF flank mode).
ORF_FLANK_FIELDS = [
    ("flank5_file_path", "5' flank  (100 bp -)",
     "FASTA of the 100 bp immediately UPSTREAM of the ATG (the '-' side, in PAM-scanning "
     "context). This lets the scan reach guide/primer positions at the start of the ORF."),
    ("flank3_file_path", "3' flank  (100 bp +)",
     "FASTA of the 100 bp immediately DOWNSTREAM of the stop codon (the '+' side, in "
     "PAM-scanning context). This lets the scan reach guide/primer positions at the end of "
     "the ORF."),
]

# The per-ORF alternative to the two flank files above: one FASTA holding
# flank5 + ORF + flank3. The ORF file is still required and locates the ORF
# inside it, so any flank length works.
ORF_PLUS_FIELD = (
    "orf_plus_file_path", "ORF + flanks",
    "FASTA holding the 5' flank, the ORF and the 3' flank as ONE sequence. The ORF "
    "file above is still required: it is found inside this sequence by exact match "
    "and the flanks are taken from either side, so the flanks may be any length. Use "
    "this when your ORF and its context live in a single file rather than three.")

# Help for the Flank inputs mode radios.
FLANK_MODE_TIPS = {
    "per_orf": (
        "Each ORF carries its own flanking sequence. Every ORF card then offers a "
        "choice: two separate flank files, or one file holding the ORF with its "
        "flanks around it. Use this when the ORFs come from different loci."),
    "global": (
        "One 5'/3' flank pair is applied to EVERY ORF, given once above the ORF "
        "cards as files or as typed sequences. Use this when all ORFs are inserted "
        "into the same locus and share the same context."),
}

# Global flank inputs (used only in global flank mode; one pair for every ORF).
# Each tuple: (key, section label, short name for buttons/messages, tooltip). The
# pamscan kwarg is '<key>_file_path' or '<key>_sequence' depending on the source.
GLOBAL_FLANK_FIELDS = [
    ("flank5", "Global 5' flank  (100 bp -)", "5' flank",
     "A single 5' flank (100 bp upstream of the ATG) applied to EVERY ORF. Use when all "
     "ORFs share the same upstream context. Open a FASTA file, or enter the sequence directly."),
    ("flank3", "Global 3' flank  (100 bp +)", "3' flank",
     "A single 3' flank (100 bp downstream of the stop) applied to EVERY ORF. Use when "
     "all ORFs share the same downstream context. Open a FASTA file, or enter the sequence directly."),
]

# Help shown on a global flank's sequence box.
SEQUENCE_TIP = (
    "Type or paste the flank sequence (A, C, G, T, or N). A pasted FASTA record, line breaks, "
    "and base numbering from a sequence viewer are all accepted -- they are stripped. The base "
    "count is shown to the right and reads 'invalid' if a non-DNA character is present.")

# Shared file inputs (apply to every ORF).
SHARED_FILE_FIELDS = [
    ("local_genome_file_path", "Genome sequence",
     "The yeast host genome (FASTA) used for off-target evaluation. PAM scanning is always "
     "performed in yeast -- the ORF is ported in from its source organism -- so this is always "
     "a yeast genome. Defaults to the bundled BY4741 genome; Browse to choose a different "
     "yeast species, strain, or variant."),
    ("codon_table_file_path", "Codon table (optional)",
     "Codon-usage table for the host genome. Leave unset to use the bundled yeast table "
     "(yeast_64_1_1_all_nuclear.cusp.txt)."),
]

# Each tuple: (kwarg key, label, default, tooltip) -- shared string settings.
STRING_FIELDS = [
    ("guidePrimerForwardSuffix", "Guide primer: forward suffix",
     "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAG",
     "Forward primer sequence that amplifies the CRISPR plasmid starting after the unique "
     "targeting sequence (UnTS). It is prepended with each 20 bp guide to build the guide "
     "plasmid via the around-the-world protocol."),
    ("insertPrimerForwardSuffix", "Insert primer: forward suffix",
     "GAAGATGTTGTCTGTTGCTCTATGTCATAT",
     "The 5' to 3' top-strand insertion-primer suffix concatenated downstream of the PAM-scan "
     "chimera-site genome homology. With the reverse suffix below it amplifies the chimeric "
     "insert for each PAM-scan site."),
    ("insertPrimerReverseSuffix", "Insert primer: reverse suffix",
     "CTTCTACAACAGACAACGAGATACAGTATA",
     "The 3' to 5' bottom-strand insertion-primer suffix concatenated upstream of the PAM-scan "
     "chimera-site genome homology. With the forward suffix above it amplifies the chimeric "
     "insert for each PAM-scan site."),
]

# Each tuple: (kwarg key, label, default, tooltip)
INT_FIELDS = [
    ("primerLength", "Primer length (bp)", 100,
     "Total length of the DNA primers used to amplify the chimera insert."),
    ("maxPamCutGap", "Max gap between two PAMs (bp)", 60,
     "Sequential PAM sites separated by more than this distance create PAM-scanning gaps with "
     "lower editing efficiency. The default of 60 reflects the empirical 30 bp rule (chimeric "
     "insertion >30 bp from the Cas9 cut site is less efficient)."),
    ("codonsSamplingGap", "Codon sampling frequency", 1,
     "Insert the chimera at every Nth codon. 1 = exhaustively after every codon, 2 = every "
     "other codon, and so on. Ignored when a codon selection file is provided."),
    ("pamInclusionThreshold", "Max PAM inclusions", 5,
     "Slightly shorter guides ('PAM inclusions') may still enable cutting; these are tracked "
     "across solutions. 5 generates ample guide designs while keeping PAM inclusions minimal."),
    ("pamInclusionSequenceThreshold", "Max PAM inclusion length (bp)", 15,
     "The minimum matched length counted as a PAM inclusion. With the default of 15, candidate "
     "inclusions longer than 15 bp are counted."),
]


# ---------------------------------------------------------------------------
# Hover tooltip
# ---------------------------------------------------------------------------

class Tooltip:
    """A lightweight hover tooltip that pops up next to a widget after a delay."""

    def __init__(self, widget, text, *, delay=350, wraplength=380):
        self.widget = widget
        self.text = text
        self.delay = delay
        self.wraplength = wraplength
        self._after_id = None
        self._tip = None
        self._anchor = widget   # the widget the box is positioned under
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, event=None):
        # Anchor the box under whichever widget the pointer entered (the input
        # field itself when hovering it), so the help reads as belonging to it.
        self._anchor = getattr(event, "widget", None) or self.widget
        self._cancel()
        self._after_id = self.widget.after(self.delay, self._show)

    def _cancel(self):
        if self._after_id is not None:
            self.widget.after_cancel(self._after_id)
            self._after_id = None

    def _show(self):
        if self._tip or not self.text:
            return
        x = self._anchor.winfo_rootx()
        y = self._anchor.winfo_rooty() + self._anchor.winfo_height() + 6
        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        self._tip.wm_geometry("+%d+%d" % (x, y))
        self._tip.attributes("-topmost", True)
        frame = tk.Frame(self._tip, background=TIP_BORDER, bd=0)
        frame.pack()
        label = tk.Label(
            frame, text=self.text, justify="left", background=TIP_BG, foreground=TEXT,
            wraplength=self.wraplength, padx=10, pady=7,
        )
        label.pack(padx=1, pady=1)

    def _hide(self, _event=None):
        self._cancel()
        if self._tip is not None:
            self._tip.destroy()
            self._tip = None


# Sequence-viewer colours (a dark "editor" panel, like the af3 residue picker).
SEQ_BG = "#0f1b28"
SEQ_FG = "#cdd9e5"
SEQ_GUTTER = "#5a6b7a"
SEQ_SELECTED = "#2f5d8a"


class CodonPicker(tk.Toplevel):
    """Modal helper to pick insertion codons graphically from the protein sequence.

    Shows the ORF's translated protein as a numbered, clickable grid (the same
    interaction as the af3 residue picker): single-click selects one residue,
    Shift/Cmd/Ctrl-click toggles residues, and positions can also be typed as
    ``"52, 89, 100-105"``. On **Use these codons** the sorted 1-based positions are
    stored in :attr:`result`; **Cancel** leaves it ``None``. Codon number *n*
    corresponds to residue *n* of the protein (codon 1 = the start Met).
    """

    def __init__(self, parent, protein, gene, initial=()):
        super().__init__(parent)
        self.protein = protein
        self.selected = {p for p in initial if 1 <= p <= len(protein)}
        self.result = None
        # Drag-select state (a left-click-drag selects a contiguous run of codons).
        self._drag_anchor = None
        self._drag_base = set()
        self._drag_last = None
        self._dragged = False

        self.title("Pick insertion codons — %s" % (gene or "ORF"))
        self.configure(bg=BG)
        self.transient(parent)
        self.resizable(True, True)
        self.minsize(640, 460)

        fam = tkfont.nametofont("TkDefaultFont").actual("family")
        avail = set(tkfont.families(self))
        mono = next((f for f in ("Menlo", "Consolas", "DejaVu Sans Mono", "Courier New")
                     if f in avail), tkfont.nametofont("TkFixedFont").actual("family"))
        self._f_label = tkfont.Font(family=fam, size=11)
        self._f_seq = tkfont.Font(family=mono, size=14)

        header = ttk.Label(
            self, style="Field.TLabel", background=BG,
            text="%s — %d residues.  Click a codon to add or remove it, or drag to select a "
                 "run (or type positions below). Picked codons are listed underneath."
                 % (gene or "ORF", len(protein)))
        header.pack(fill="x", padx=12, pady=(12, 6))

        viewer = ttk.Frame(self, style="TFrame")
        viewer.pack(fill="both", expand=True, padx=12)
        self.view = tk.Text(viewer, wrap="none", font=self._f_seq, height=12,
                            background=SEQ_BG, foreground=SEQ_FG, insertwidth=0,
                            cursor="arrow", spacing1=2, spacing3=2, borderwidth=0,
                            highlightthickness=0)
        vbar = ttk.Scrollbar(viewer, orient="vertical", command=self.view.yview)
        self.view.configure(yscrollcommand=vbar.set)
        self.view.grid(row=0, column=0, sticky="nsew")
        vbar.grid(row=0, column=1, sticky="ns")
        viewer.rowconfigure(0, weight=1)
        viewer.columnconfigure(0, weight=1)
        # Read-only but clickable: swallow typing/paste, keep tag bindings live.
        for evt in ("<Key>", "<<Paste>>", "<<Cut>>", "<<Clear>>"):
            self.view.bind(evt, lambda e: "break")
        self.view.tag_configure("gutter", foreground=SEQ_GUTTER)
        self.view.tag_configure("selected", background=SEQ_SELECTED)
        # Press on a residue starts a click (toggle) or a drag (run select); the
        # motion/release are handled at the widget level so a drag can sweep across
        # the whole grid, not just the residue first pressed.
        for evt in ("<Button-1>", "<Shift-Button-1>", "<Command-Button-1>", "<Control-Button-1>"):
            self.view.tag_bind("residue", evt, self._on_press)
        self.view.bind("<B1-Motion>", self._on_drag)
        self.view.bind("<ButtonRelease-1>", self._on_release)
        self.view.tag_bind("residue", "<Enter>",
                           lambda e: self.view.configure(cursor="hand2"))
        self.view.tag_bind("residue", "<Leave>",
                           lambda e: self.view.configure(cursor="arrow"))
        # Wheel/trackpad over the grid scrolls the grid itself (and, via the "break",
        # never the window behind the picker).
        self._wheel_bind(self.view, self.view)

        ctl = ttk.Frame(self, style="TFrame")
        ctl.pack(fill="x", padx=12, pady=(6, 2))
        ttk.Label(ctl, style="Field.TLabel", background=BG,
                  text="Type positions:").pack(side="left")
        self.pos_entry = ttk.Entry(ctl, width=20)
        self.pos_entry.pack(side="left", padx=4)
        self.pos_entry.bind("<Return>", lambda e: self._add_typed())
        ttk.Button(ctl, text="Add", style="Small.TButton", width=6,
                   command=self._add_typed).pack(side="left")
        ttk.Button(ctl, text="Select all", style="Small.TButton",
                   command=self._select_all).pack(side="left", padx=4)
        ttk.Button(ctl, text="Clear", style="Small.TButton", width=6,
                   command=self._clear).pack(side="left")

        # Picked-codons list: one line per contiguous run, each individually removable.
        list_head = ttk.Frame(self, style="TFrame")
        list_head.pack(fill="x", padx=12, pady=(8, 0))
        self.sel_label = ttk.Label(list_head, style="Field.TLabel", background=BG,
                                   foreground=HEADER_BG, text="Picked codons: none")
        self.sel_label.pack(side="left")

        list_wrap = ttk.Frame(self, style="TFrame")
        list_wrap.pack(fill="x", padx=12, pady=(2, 2))
        self.list_canvas = tk.Canvas(list_wrap, height=140, background=CARD,
                                     highlightthickness=1, highlightbackground="#c9d3de")
        list_sb = ttk.Scrollbar(list_wrap, orient="vertical", command=self.list_canvas.yview)
        self.list_canvas.configure(yscrollcommand=list_sb.set)
        self.list_canvas.pack(side="left", fill="both", expand=True)
        list_sb.pack(side="right", fill="y")
        self.list_inner = ttk.Frame(self.list_canvas, style="Card.TFrame")
        self._list_window = self.list_canvas.create_window((0, 0), window=self.list_inner, anchor="nw")
        self.list_inner.bind(
            "<Configure>",
            lambda e: self.list_canvas.configure(scrollregion=self.list_canvas.bbox("all")))
        self.list_canvas.bind(
            "<Configure>",
            lambda e: self.list_canvas.itemconfigure(self._list_window, width=e.width))
        # Wheel/trackpad scrolling, bound per-widget (never bind_all — the main GUI
        # relies on its own bind_all scroll handlers).
        self._wheel_bind(self.list_canvas, self.list_canvas)
        self._wheel_bind(self.list_inner, self.list_canvas)

        actions = ttk.Frame(self, style="TFrame")
        actions.pack(fill="x", padx=12, pady=(4, 12))
        ttk.Button(actions, text="Use these codons", style="Browse.TButton",
                   command=self._accept).pack(side="right")
        ttk.Button(actions, text="Cancel", style="Small.TButton",
                   command=self._cancel).pack(side="right", padx=8)

        self._render()
        self.bind("<Escape>", lambda e: self._cancel())
        self.protocol("WM_DELETE_WINDOW", self._cancel)
        self.update_idletasks()
        self.grab_set()
        self.pos_entry.focus_set()

    # --- rendering -----------------------------------------------------------
    def _render(self):
        """Draw the protein as a position-ruled grid; each residue is its own tag."""
        txt = self.view
        txt.delete("1.0", "end")
        seq, per_line, group = self.protein, 30, 10
        for start in range(0, len(seq), per_line):
            chunk = seq[start:start + per_line]
            txt.insert("end", "%6d  " % (start + 1), ("gutter",))
            for k, ch in enumerate(chunk):
                pos = start + k + 1
                tags = ("residue", "res_%d" % pos)
                txt.insert("end", ch, tags)   # residue letter + trailing space share the
                txt.insert("end", " ", tags)  # tag: a wider, easier click target per residue
                if (k + 1) % group == 0 and (k + 1) < len(chunk):
                    txt.insert("end", "  ", ("gutter",))
            txt.insert("end", " %d\n" % (start + len(chunk)), ("gutter",))
        self._changed()

    def _residue_at(self, event):
        idx = self.view.index("@%d,%d" % (event.x, event.y))
        for probe in (idx, "%s-1c" % idx):
            for t in self.view.tag_names(probe):
                if t.startswith("res_"):
                    return int(t[4:])
        return None

    def _on_press(self, event):
        """Press on a residue: remember the anchor; a plain click toggles it (on
        release), a drag from here selects the run it sweeps."""
        self._drag_anchor = self._residue_at(event)
        self._drag_base = set(self.selected)
        self._drag_last = self._drag_anchor
        self._dragged = False
        return "break"

    def _on_drag(self, event):
        """Left-click-drag: add every codon between the anchor and the pointer."""
        if self._drag_anchor is None:
            return "break"
        pos = self._residue_at(event)
        if pos is None:                 # over a gap/gutter: hold the last residue
            pos = self._drag_last
        if pos is None:
            return "break"
        self._drag_last = pos
        if pos != self._drag_anchor:
            self._dragged = True
        lo, hi = sorted((self._drag_anchor, pos))
        self.selected = self._drag_base | set(range(lo, hi + 1))
        # Keep the sweep responsive: update the grid + count now, defer the (heavier)
        # list rebuild to release.
        self._refresh_highlight()
        self._update_label()
        return "break"

    def _on_release(self, event):
        if self._drag_anchor is None:
            return "break"
        if not self._dragged:           # no motion => a click that toggles one codon
            self.selected = set(self._drag_base) ^ {self._drag_anchor}
        self._drag_anchor = None
        self._changed()
        return "break"

    def _changed(self):
        """Single point that re-syncs the grid highlight, the count, and the list."""
        self._refresh_highlight()
        self._update_label()
        self._render_list()

    def _refresh_highlight(self):
        self.view.tag_remove("selected", "1.0", "end")
        for pos in self.selected:
            r = self.view.tag_ranges("res_%d" % pos)
            if r:
                self.view.tag_add("selected", r[0], r[1])

    def _update_label(self):
        n = len(self.selected)
        self.sel_label.configure(
            text="Picked codons: none" if not n else "Picked codons (%d)" % n)

    def _runs(self):
        """Collapse the selected positions into sorted, contiguous (lo, hi) runs."""
        runs = []
        for p in sorted(self.selected):
            if runs and p == runs[-1][1] + 1:
                runs[-1][1] = p
            else:
                runs.append([p, p])
        return [(lo, hi) for lo, hi in runs]

    def _render_list(self):
        """Rebuild the picked-codons list: one row per run, each with a Remove button."""
        for w in self.list_inner.winfo_children():
            w.destroy()
        self.list_inner.columnconfigure(0, weight=1)
        runs = self._runs()
        if not runs:
            empty = ttk.Label(self.list_inner, style="Field.TLabel", foreground=MUTED,
                              text="No codons picked yet — click residues above or type positions.")
            empty.grid(row=0, column=0, sticky="w", padx=8, pady=6)
            self._wheel_bind(empty, self.list_canvas)
            return
        for i, (lo, hi) in enumerate(runs):
            if lo == hi:
                text = "Codon %d  (%s%d)" % (lo, self.protein[lo - 1], lo)
            else:
                text = "Codons %d–%d  (%d codons)" % (lo, hi, hi - lo + 1)
            lbl = ttk.Label(self.list_inner, text=text, style="Field.TLabel")
            lbl.grid(row=i, column=0, sticky="w", padx=(8, 6), pady=1)
            btn = ttk.Button(self.list_inner, text="Remove", style="Small.TButton", width=8,
                             command=lambda a=lo, b=hi: self._remove_run(a, b))
            btn.grid(row=i, column=1, sticky="e", padx=(6, 8), pady=1)
            self._wheel_bind(lbl, self.list_canvas)

    def _remove_run(self, lo, hi):
        self.selected -= set(range(lo, hi + 1))
        self._changed()

    # --- scrolling ------------------------------------------------------------
    # Bindings are per-widget (never bind_all) and return "break", so a scroll inside
    # the picker scrolls its own widget and never the main window behind it. Direction
    # matches the main window (reversed from the raw platform sign).
    def _wheel_bind(self, widget, target):
        """Scroll *target* when the wheel/trackpad is used over *widget*."""
        widget.bind("<MouseWheel>", lambda e: self._scroll_wheel(target, e))
        widget.bind("<Button-4>", lambda e: self._scroll_wheel(target, e))
        widget.bind("<Button-5>", lambda e: self._scroll_wheel(target, e))
        try:
            widget.bind("<TouchpadScroll>", lambda e: self._scroll_touch(target, e))
        except tk.TclError:
            pass

    def _scroll_wheel(self, target, event):
        if getattr(event, "num", None) == 4:
            step = 1
        elif getattr(event, "num", None) == 5:
            step = -1
        else:
            step = 1 if getattr(event, "delta", 0) > 0 else -1
        target.yview_scroll(step, "units")
        return "break"

    def _scroll_touch(self, target, event):
        try:
            _dx, dy = self.tk.call("tk::PreciseScrollDeltas", event.delta)
            if dy:
                target.yview_scroll(int(dy), "units")
        except tk.TclError:
            pass
        return "break"

    # --- controls ------------------------------------------------------------
    def _add_typed(self):
        added = parse_codon_positions(self.pos_entry.get(), len(self.protein))
        if not added:
            messagebox.showerror("No positions",
                                 "Enter positions like '52, 89, 100-105'.", parent=self)
            return
        self.selected |= set(added)
        self.pos_entry.delete(0, "end")
        self._changed()

    def _select_all(self):
        self.selected = set(range(1, len(self.protein) + 1))
        self._changed()

    def _clear(self):
        self.selected = set()
        self._changed()

    def _accept(self):
        self.result = sorted(self.selected)
        self.grab_release()
        self.destroy()

    def _cancel(self):
        self.result = None
        self.grab_release()
        self.destroy()


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

def main():
    root = tk.Tk()
    root.title(f"PAM Scanning v{__version__}")
    root.configure(bg=BG)
    root.minsize(1380, 860)
    root.geometry("1580x920")

    # Fonts.
    base = tkfont.nametofont("TkDefaultFont")
    family = base.actual("family")
    # Pick a monospace family that exists on this platform; fall back to Tk's
    # built-in fixed font (which Tk maps to a sensible default everywhere).
    available = set(tkfont.families(root))
    mono_family = next(
        (f for f in ("Menlo", "Consolas", "DejaVu Sans Mono", "Courier New")
         if f in available),
        tkfont.nametofont("TkFixedFont").actual("family"),
    )
    f_title = tkfont.Font(family=family, size=20, weight="bold")
    f_subtitle = tkfont.Font(family=family, size=11)
    f_section = tkfont.Font(family=family, size=12, weight="bold")
    f_label = tkfont.Font(family=family, size=11)
    f_mono = tkfont.Font(family=mono_family, size=10)
    f_button = tkfont.Font(family=family, size=11)
    f_run = tkfont.Font(family=family, size=13, weight="bold")

    # Theme + styles.
    style = ttk.Style()
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass
    style.configure("TFrame", background=BG)
    style.configure("Card.TFrame", background=CARD)
    style.configure("Header.TFrame", background=HEADER_BG)
    style.configure("Title.TLabel", background=HEADER_BG, foreground=HEADER_FG, font=f_title)
    style.configure("Subtitle.TLabel", background=HEADER_BG, foreground=SUBTITLE_FG, font=f_subtitle)
    style.configure("Section.TLabel", background=BG, foreground=HEADER_BG, font=f_section)
    style.configure("Field.TLabel", background=CARD, foreground=TEXT, font=f_label)
    style.configure("OrfHead.TLabel", background=CARD, foreground=HEADER_BG, font=f_section)
    style.configure("Path.TLabel", background=CARD, foreground=MUTED, font=f_label)
    style.configure("Status.TLabel", background=BG, foreground=MUTED, font=f_subtitle)
    style.configure("Console.TFrame", background=CONSOLE_BG)
    style.configure("ConsoleHead.TFrame", background=HEADER_BG)
    style.configure("ConsoleHead.TLabel", background=HEADER_BG, foreground=HEADER_FG, font=f_section)
    style.configure("Mode.TRadiobutton", background=CARD, foreground=TEXT, font=f_label)
    style.map("Mode.TRadiobutton", background=[("active", CARD)])
    style.configure("TEntry", fieldbackground="#ffffff", padding=4)
    style.configure("Browse.TButton", font=f_button, padding=(10, 4))
    style.configure("Small.TButton", font=f_button, padding=(8, 3))
    style.configure(
        "Accent.TButton", font=f_run, foreground="#ffffff", background=ACCENT,
        padding=(26, 12), borderwidth=0,
    )
    style.map(
        "Accent.TButton",
        background=[("active", ACCENT_ACTIVE), ("disabled", "#9bbf9d")],
        foreground=[("disabled", "#eef2f5")],
    )

    # Integer-entry validation.
    def _validate_int(proposed):
        return proposed == "" or proposed.isdigit()

    vint = (root.register(_validate_int), "%P")

    # Shared state variables.
    shared_file_vars = {key: tk.StringVar(value=PLACEHOLDER) for key, _, _ in SHARED_FILE_FIELDS}
    # Default the genome to the bundled yeast genome (decompressed on first use).
    try:
        from pam_scanning.chimeras import default_genome_path
        shared_file_vars["local_genome_file_path"].set(default_genome_path())
    except Exception:   # never let a packaging hiccup stop the GUI from opening
        pass
    # Each global flank has a source ("file" or "sequence") and a variable per source.
    global_flank_source = {key: tk.StringVar(value="file") for key, _, _, _ in GLOBAL_FLANK_FIELDS}
    global_flank_path_vars = {key: tk.StringVar(value=PLACEHOLDER) for key, _, _, _ in GLOBAL_FLANK_FIELDS}
    global_flank_seq_vars = {key: tk.StringVar(value="") for key, _, _, _ in GLOBAL_FLANK_FIELDS}
    global_flank_widgets = {}                    # key -> {"file": [...], "sequence": [...]}
    AUTO_DB = "Auto — built from the genome on first use"
    blast_db_var = tk.StringVar(value=AUTO_DB)   # AUTO_DB => build from the genome; else -db override
    string_vars = {key: tk.StringVar(value=default) for key, _, default, _ in STRING_FIELDS}
    int_vars = {key: tk.StringVar(value=str(default)) for key, _, default, _ in INT_FIELDS}
    output_var = tk.StringVar(value=os.getcwd())
    status_var = tk.StringVar(value="Ready.")
    flank_mode = tk.StringVar(value="per_orf")   # "per_orf" or "global"
    orf_entries = []   # one dict of StringVars (+ widgets) per queued ORF

    # Every file/folder dialog opens at the last directory used, remembered across
    # sessions. A single shared memory means the app stays wherever you last were.
    dialog_dir = {"path": _load_last_dir()}

    # Progress console: worker threads enqueue text; the UI thread drains it.
    console_queue = queue.Queue()

    def console_log(text):
        """Thread-safe: queue *text* for the progress console (no trailing newline added)."""
        console_queue.put(text)

    def console_banner(text):
        """Queue a highlighted banner line (e.g. a stage or ORF header)."""
        console_queue.put((text, "banner"))

    def console_error(text):
        """Queue an error-styled line."""
        console_queue.put((text, "error"))

    def remember_dir(chosen):
        """After a dialog picks `chosen`, point the next dialog at its directory."""
        if not chosen:
            return
        directory = chosen if os.path.isdir(chosen) else os.path.dirname(chosen)
        if directory:
            dialog_dir["path"] = directory
            _save_last_dir(directory)

    # --- Scrollable body so the form fits any screen ---------------------
    outer = ttk.Frame(root, style="TFrame")
    outer.pack(fill="both", expand=True)

    # Split the window: the scrollable form on the left, a live console on the right.
    paned = ttk.PanedWindow(outer, orient="horizontal")
    paned.pack(fill="both", expand=True)

    form_pane = ttk.Frame(paned, style="TFrame")
    canvas = tk.Canvas(form_pane, background=BG, highlightthickness=0, yscrollincrement=1)
    scrollbar = ttk.Scrollbar(form_pane, orient="vertical", command=canvas.yview)
    body = ttk.Frame(canvas, style="TFrame")
    body.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
    body_window = canvas.create_window((0, 0), window=body, anchor="nw")
    canvas.bind("<Configure>", lambda e: canvas.itemconfigure(body_window, width=e.width))
    canvas.configure(yscrollcommand=scrollbar.set)
    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")
    paned.add(form_pane, weight=3)

    # Scrolling. Tk 9 delivers trackpad two-finger scroll as <TouchpadScroll> (NOT
    # <MouseWheel>); its precise deltas are decoded with tk::PreciseScrollDeltas. With
    # yscrollincrement=1 a canvas "unit" is one pixel, so scrolling is smooth.
    def _touchpad_dy(event):
        try:
            _dx, dy = (int(v) for v in root.tk.call("tk::PreciseScrollDeltas", event.delta))
            return -dy   # reversed so the main window scrolls the way the trackpad expects
        except Exception:
            return 0

    def _wheel_direction(event):
        """Mouse-wheel direction as +1 (content down) / -1 (content up) / 0. Reversed
        from the raw platform sign so the main window scrolls the same way as the
        codon picker."""
        if event.num == 4:
            return 1
        if event.num == 5:
            return -1
        if not event.delta:
            return 0
        return -1 if event.delta < 0 else 1

    def _event_in_main(event):
        """True only when the scroll happened over the main window, not a modal on top
        of it (e.g. the codon picker). The form's scroll handlers are bound with
        ``bind_all``, so without this guard scrolling inside the picker would also move
        the form underneath it."""
        w = getattr(event, "widget", None)
        try:
            return w is not None and w.winfo_toplevel() is root
        except (AttributeError, tk.TclError):
            return False

    def _scroll_pixels(widget, accum, key, pixels, what):
        """Scroll *widget* by a possibly-fractional pixel amount, banking the remainder."""
        accum[key] += pixels
        steps = int(accum[key])
        if steps:
            accum[key] -= steps
            widget.yview_scroll(steps, what)

    form_accum = {"y": 0.0}

    def _on_touchpad(event):
        if not _event_in_main(event):
            return
        dy = _touchpad_dy(event)
        if dy:
            _scroll_pixels(canvas, form_accum, "y", dy * TOUCHPAD_SCALE, "units")

    def _on_mousewheel(event):
        if not _event_in_main(event):
            return
        direction = _wheel_direction(event)
        if direction:
            canvas.yview_scroll(direction * MOUSE_WHEEL_PIXELS, "units")

    try:
        canvas.bind_all("<TouchpadScroll>", _on_touchpad)   # Tk >= 8.7 / 9
    except tk.TclError:
        pass   # older Tk: the trackpad arrives as <MouseWheel>, handled below
    for _seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
        canvas.bind_all(_seq, _on_mousewheel)

    # --- Progress console (right pane) -----------------------------------
    console_pane = ttk.Frame(paned, style="Console.TFrame")
    console_head = ttk.Frame(console_pane, style="ConsoleHead.TFrame")
    console_head.pack(fill="x")
    ttk.Label(console_head, text="Progress", style="ConsoleHead.TLabel").pack(
        side="left", padx=14, pady=9)

    def clear_console():
        console.configure(state="normal")
        console.delete("1.0", "end")
        console.configure(state="disabled")

    ttk.Button(console_head, text="Clear", style="Small.TButton", command=clear_console).pack(
        side="right", padx=10, pady=7)

    console_wrap = ttk.Frame(console_pane, style="Console.TFrame")
    console_wrap.pack(fill="both", expand=True)
    console = tk.Text(
        console_wrap, background=CONSOLE_BG, foreground=CONSOLE_FG, insertbackground=CONSOLE_FG,
        font=f_mono, wrap="word", relief="flat", borderwidth=0, padx=12, pady=10,
        width=48, height=20, state="disabled", highlightthickness=0,
    )
    console_scroll = ttk.Scrollbar(console_wrap, orient="vertical", command=console.yview)
    console.configure(yscrollcommand=console_scroll.set)
    console_scroll.pack(side="right", fill="y")
    console.pack(side="left", fill="both", expand=True)
    console.tag_configure("banner", foreground="#8fd3ff")
    console.tag_configure("muted", foreground=CONSOLE_MUTED)
    console.tag_configure("error", foreground="#ff9a8f")
    paned.add(console_pane, weight=2)

    console_images = []   # keep PhotoImage refs so Tk doesn't garbage-collect them

    def embed_plot(png_path):
        """Insert a saved plot PNG into the console (main thread), scaled to the pane."""
        try:
            img = tk.PhotoImage(file=png_path)
        except tk.TclError:
            return   # e.g. Tk built without PNG support; the file is still saved to QC
        pane_width = max(360, console.winfo_width() - 28)
        factor = max(1, -(-img.width() // pane_width))   # ceil division
        if factor > 1:
            img = img.subsample(factor, factor)
        console_images.append(img)
        console.configure(state="normal")
        console.insert("end", "\n")
        console.image_create("end", image=img)
        console.insert("end", "\n")
        console.see("end")
        console.configure(state="disabled")

    def _append_console(text, tag=None):
        console.configure(state="normal")
        console.insert("end", text, tag or ())
        console.see("end")
        console.configure(state="disabled")

    # Scroll the console over itself; return "break" so it doesn't also move the form.
    console_accum = {"y": 0.0}

    def _console_touchpad(event):
        dy = _touchpad_dy(event)
        if dy:
            _scroll_pixels(console, console_accum, "y", dy * TOUCHPAD_SCALE, "pixels")
        return "break"

    def _console_wheel(event):
        direction = _wheel_direction(event)
        if direction:
            console.yview_scroll(direction * MOUSE_WHEEL_PIXELS, "pixels")
        return "break"

    try:
        console.bind("<TouchpadScroll>", _console_touchpad)   # Tk >= 8.7 / 9
    except tk.TclError:
        pass
    console.bind("<MouseWheel>", _console_wheel)
    console.bind("<Button-4>", _console_wheel)
    console.bind("<Button-5>", _console_wheel)

    def drain_console():
        try:
            while True:
                item = console_queue.get_nowait()
                if isinstance(item, tuple):
                    _append_console(item[0], item[1])
                else:
                    _append_console(item)
        except queue.Empty:
            pass
        root.after(120, drain_console)

    drain_console()

    # --- Header banner ---------------------------------------------------
    header = ttk.Frame(body, style="Header.TFrame")
    header.pack(fill="x")
    inner_h = ttk.Frame(header, style="Header.TFrame")
    inner_h.pack(fill="x", padx=28, pady=(20, 18))
    ttk.Label(inner_h, text="PAM Scanning", style="Title.TLabel").pack(anchor="w")
    ttk.Label(
        inner_h,
        text="Design CRISPR guides and primers for ORF-wide chimera insertion.",
        style="Subtitle.TLabel",
    ).pack(anchor="w", pady=(4, 0))

    content = ttk.Frame(body, style="TFrame")
    content.pack(fill="both", expand=True, padx=24, pady=18)

    def attach_tip(anchor, text, *widgets):
        """Reveal `text` as a hover tooltip over the anchor and any extra row widgets."""
        tip = Tooltip(anchor, text)
        for w in widgets:
            w.bind("<Enter>", tip._schedule, add="+")
            w.bind("<Leave>", tip._hide, add="+")
            w.bind("<ButtonPress>", tip._hide, add="+")
        return tip

    def section(title, weighted_col):
        ttk.Label(content, text=title, style="Section.TLabel").pack(anchor="w", pady=(14, 6))
        card = ttk.Frame(content, style="Card.TFrame")
        card.pack(fill="x")
        card.columnconfigure(weighted_col, weight=1)
        return card

    def add_file_row(card, row, button_label, tip, var, *, button_text=None, button_width=26,
                     optional=False):
        """Add a path-label + Browse-button row; return its widgets (for show/hide).

        By default the button reads "Browse  <button_label>"; pass *button_text* to
        set the label verbatim (and *button_width* to size it).

        With *optional*, the row's label reads "No file selected  (optional)" while
        nothing is chosen, so the reader does not have to hover to learn the input
        can be left out. The path variable itself is untouched and still holds the
        PLACEHOLDER sentinel that the pipeline treats as "not provided".
        """
        display = var
        if optional:
            display = tk.StringVar()

            def _mirror(*_, source=var, shown=display):
                value = source.get()
                shown.set("%s  (optional)" % PLACEHOLDER
                          if not value or value == PLACEHOLDER else value)

            _mirror()
            var.trace_add("write", _mirror)

        path_lbl = ttk.Label(card, textvariable=display, style="Path.TLabel",
                             wraplength=560, anchor="w", justify="left")
        path_lbl.grid(row=row, column=0, sticky="we", padx=(14, 6), pady=10)

        def browse(v=var):
            # Open where this field already points (e.g. the bundled genome), else
            # fall back to the shared last-used directory.
            current = v.get()
            start = (os.path.dirname(current)
                     if current and current != PLACEHOLDER and os.path.isfile(current)
                     else dialog_dir["path"])
            chosen = filedialog.askopenfilename(initialdir=start)
            if chosen:
                remember_dir(chosen)
                v.set(chosen)

        btn = ttk.Button(card, text=button_text or ("Browse  " + button_label),
                         style="Browse.TButton", command=browse, width=button_width)
        btn.grid(row=row, column=1, sticky="e", padx=(6, 14), pady=10)
        attach_tip(path_lbl, tip, btn)
        return [path_lbl, btn]

    def add_entry_row(card, row, label, tip, var, *, mono=False, validate=False):
        lbl = ttk.Label(card, text=label, style="Field.TLabel")
        lbl.grid(row=row, column=0, sticky="w", padx=(14, 6), pady=8)
        kwargs = {"textvariable": var, "font": f_mono if mono else f_label}
        if validate:
            kwargs["validate"] = "key"
            kwargs["validatecommand"] = vint
            kwargs["justify"] = "right"
            kwargs["width"] = 12
        entry = ttk.Entry(card, **kwargs)
        if validate:
            entry.grid(row=row, column=1, sticky="e", padx=(6, 14), pady=8)
        else:
            entry.grid(row=row, column=1, sticky="we", padx=(6, 14), pady=8)
        attach_tip(lbl, tip, entry)
        return entry

    # --- Flank inputs mode ----------------------------------------------
    ttk.Label(content, text="Flank inputs", style="Section.TLabel").pack(anchor="w", pady=(14, 6))
    mode_card = ttk.Frame(content, style="Card.TFrame")
    mode_card.pack(fill="x")
    radios = ttk.Frame(mode_card, style="Card.TFrame")
    radios.pack(anchor="w", padx=14, pady=10)
    per_orf_radio = ttk.Radiobutton(radios, text="Per-ORF flanks", value="per_orf",
                                    variable=flank_mode, style="Mode.TRadiobutton")
    per_orf_radio.pack(side="left", padx=(0, 24))
    global_radio = ttk.Radiobutton(radios, text="Global flanks (one 5'/3' pair for all ORFs)",
                                   value="global", variable=flank_mode, style="Mode.TRadiobutton")
    global_radio.pack(side="left")
    attach_tip(per_orf_radio, FLANK_MODE_TIPS["per_orf"])
    attach_tip(global_radio, FLANK_MODE_TIPS["global"])

    # Global flank pickers (shown only in global mode). The card stays packed
    # inside the holder for the life of the window and the HOLDER is what gets
    # hidden: Tk never shrinks a frame back to nothing once it has managed a
    # slave, so forgetting the card would leave the holder standing at its old
    # height and open a blank gap above the ORFs section.
    global_flank_holder = ttk.Frame(content, style="TFrame")
    global_flank_holder.pack(fill="x")
    global_flank_card = ttk.Frame(global_flank_holder, style="Card.TFrame")
    global_flank_card.columnconfigure(0, weight=1)
    global_flank_card.pack(fill="x")

    def refresh_global_flank_source(*_):
        """Show the Browse row or the sequence box, per flank, per selected source."""
        for key, widgets in global_flank_widgets.items():
            from_file = global_flank_source[key].get() == "file"
            shown, hidden = ("file", "sequence") if from_file else ("sequence", "file")
            for w in widgets[shown]:
                w.grid()
            for w in widgets[hidden]:
                w.grid_remove()

    def add_global_flank(card, base_row, key, label, short, tip):
        """Build one global flank: a File/Sequence toggle over a Browse row and a sequence box."""
        head = ttk.Frame(card, style="Card.TFrame")
        head.grid(row=base_row, column=0, columnspan=2, sticky="we", padx=14, pady=(10, 0))
        head_lbl = ttk.Label(head, text=label, style="Field.TLabel")
        head_lbl.pack(side="left")
        source = global_flank_source[key]
        ttk.Radiobutton(head, text="Enter sequence", value="sequence", variable=source,
                        style="Mode.TRadiobutton").pack(side="right")
        ttk.Radiobutton(head, text="From file", value="file", variable=source,
                        style="Mode.TRadiobutton").pack(side="right", padx=(0, 16))
        attach_tip(head_lbl, tip)

        file_widgets = add_file_row(card, base_row + 1, short, tip, global_flank_path_vars[key])

        seq_var = global_flank_seq_vars[key]
        seq_entry = ttk.Entry(card, textvariable=seq_var, font=f_mono)
        seq_entry.grid(row=base_row + 2, column=0, sticky="we", padx=(14, 6), pady=10)
        length_lbl = ttk.Label(card, text="0 bp", style="Path.TLabel", anchor="e", width=12)
        length_lbl.grid(row=base_row + 2, column=1, sticky="e", padx=(6, 14), pady=10)
        attach_tip(seq_entry, SEQUENCE_TIP, length_lbl)

        def update_length(*_):
            """Echo the parsed base count so a mistyped/short flank is obvious before running."""
            text = seq_var.get().strip()
            if not text:
                length_lbl.config(text="0 bp", foreground=MUTED)
                return
            try:
                length_lbl.config(text="%d bp" % len(parse_sequence_text(text)), foreground=MUTED)
            except ValueError:
                length_lbl.config(text="invalid", foreground=INVALID)

        seq_var.trace_add("write", update_length)
        source.trace_add("write", refresh_global_flank_source)
        global_flank_widgets[key] = {"file": file_widgets, "sequence": [seq_entry, length_lbl]}

    for i, (key, label, short, tip) in enumerate(GLOBAL_FLANK_FIELDS):
        add_global_flank(global_flank_card, i * 3, key, label, short, tip)
    refresh_global_flank_source()

    # --- ORFs (one or more) ---------------------------------------------
    ttk.Label(content, text="ORFs", style="Section.TLabel").pack(anchor="w", pady=(14, 6))
    orf_container = ttk.Frame(content, style="TFrame")
    orf_container.pack(fill="x")

    def renumber_orfs():
        for i, entry in enumerate(orf_entries, start=1):
            entry["head"].config(text="ORF %d" % i)
            entry["remove_btn"].config(state=("disabled" if len(orf_entries) == 1 else "normal"))

    def remove_orf(entry):
        if len(orf_entries) == 1:
            return
        entry["card"].destroy()
        orf_entries.remove(entry)
        renumber_orfs()

    def apply_orf_flank_source(entry):
        """Within an ORF card, show either the two flank rows or the combined row."""
        if flank_mode.get() != "per_orf":
            return
        separate = entry["flank_source"].get() == "separate"
        for w in entry["flank_widgets"]:
            w.grid() if separate else w.grid_remove()
        for w in entry["orf_plus_widgets"]:
            w.grid_remove() if separate else w.grid()

    def apply_flank_mode_to_entry(entry):
        """Global mode hides everything per-ORF about flanks, including the choice."""
        per_orf = flank_mode.get() == "per_orf"
        if not per_orf:
            for w in (entry["source_widgets"] + entry["flank_widgets"]
                      + entry["orf_plus_widgets"]):
                w.grid_remove()
            return
        for w in entry["source_widgets"]:
            w.grid()
        apply_orf_flank_source(entry)

    def update_picked_label(entry):
        positions = entry.get("codon_positions") or []
        var = entry["_picked_var"]
        if not positions:
            var.set(NO_CODONS_PICKED)
            return
        shown = ", ".join(str(p) for p in positions[:14]) + (" …" if len(positions) > 14 else "")
        var.set("%d codon(s) picked: %s" % (len(positions), shown))

    def open_codon_picker(entry):
        orf_path = entry["orf_file_path"].get()
        if not _is_set(orf_path) or not os.path.isfile(orf_path):
            messagebox.showerror(
                "Choose an ORF first",
                "Select this ORF's FASTA file before picking codons — the protein "
                "sequence is translated from it.")
            return
        try:
            from pam_scanning.chimeras import _read_fasta_sequence
            dna = _read_fasta_sequence(orf_path)
            protein = fetch_cds.translate(dna)
        except Exception as exc:
            messagebox.showerror("Could not read ORF", "Failed to read the ORF FASTA:\n%s" % exc)
            return
        if not protein:
            messagebox.showerror(
                "Empty protein",
                "The ORF did not translate to any residues. Check that it is a coding "
                "DNA sequence (ATG…stop).")
            return
        gene = entry["geneName"].get().strip() or gene_name_from_orf_path(orf_path)
        picker = CodonPicker(root, protein, gene, entry.get("codon_positions") or [])
        root.wait_window(picker)
        if picker.result is not None:      # None = cancelled; [] = deliberately cleared
            entry["codon_positions"] = picker.result
            update_picked_label(entry)

    def add_orf():
        card = ttk.Frame(orf_container, style="Card.TFrame")
        card.pack(fill="x", pady=(0, 10))
        card.columnconfigure(0, weight=1)

        head_bar = ttk.Frame(card, style="Card.TFrame")
        head_bar.grid(row=0, column=0, columnspan=2, sticky="we", padx=14, pady=(8, 0))
        head = ttk.Label(head_bar, text="ORF", style="OrfHead.TLabel")
        head.pack(side="left")
        entry = {"card": card, "head": head}
        remove_btn = ttk.Button(head_bar, text="Remove", style="Small.TButton",
                                command=lambda e=entry: remove_orf(e))
        remove_btn.pack(side="right")
        entry["remove_btn"] = remove_btn

        gene_var = tk.StringVar(value="")
        entry["geneName"] = gene_var
        add_entry_row(card, 1, "Gene name", "A short gene-name label used in this ORF's output "
                      "file names.", gene_var)

        orf_var = tk.StringVar(value=PLACEHOLDER)
        entry["orf_file_path"] = orf_var
        add_file_row(card, 2, "ORF", "Open the open reading frame (ORF) FASTA file for this "
                     "gene. The ORF should begin with the ATG start codon and end with a stop "
                     "codon. If the gene name is left blank it is derived from this file name.",
                     orf_var)

        def _autofill_gene(*_, g=gene_var, o=orf_var):
            # When an ORF file is chosen and no gene name is set, derive one from it.
            if not g.get().strip():
                path = o.get()
                if path and path != PLACEHOLDER and os.path.isfile(path):
                    g.set(gene_name_from_orf_path(path))

        orf_var.trace_add("write", _autofill_gene)

        # How this ORF's flanks arrive: two separate files, or one file holding
        # the ORF with its flanks around it. Row 3; hidden in global flank mode.
        source_var = tk.StringVar(value="separate")
        entry["flank_source"] = source_var
        source_bar = ttk.Frame(card, style="Card.TFrame")
        source_bar.grid(row=3, column=0, columnspan=2, sticky="we", padx=14, pady=(4, 0))
        ttk.Label(source_bar, text="Flanks:", style="Path.TLabel").pack(side="left", padx=(0, 10))
        sep_radio = ttk.Radiobutton(source_bar, text="Separate flank files", value="separate",
                                    variable=source_var, style="Mode.TRadiobutton")
        sep_radio.pack(side="left", padx=(0, 18))
        plus_radio = ttk.Radiobutton(source_bar, text="ORF + flanks in one file", value="combined",
                                     variable=source_var, style="Mode.TRadiobutton")
        plus_radio.pack(side="left")
        attach_tip(sep_radio, "Give this ORF's 5' and 3' flanks as two separate FASTA files.")
        attach_tip(plus_radio, ORF_PLUS_FIELD[2])
        entry["source_widgets"] = [source_bar]

        # Separate flank rows (4-5).
        flank_widgets = []
        for i, (key, blabel, tip) in enumerate(ORF_FLANK_FIELDS):
            var = tk.StringVar(value=PLACEHOLDER)
            entry[key] = var
            flank_widgets += add_file_row(card, 4 + i, blabel, tip, var)
        entry["flank_widgets"] = flank_widgets

        # Combined ORF + flanks row (6); the alternative to the two rows above.
        plus_key, plus_label, plus_tip = ORF_PLUS_FIELD
        plus_var = tk.StringVar(value=PLACEHOLDER)
        entry[plus_key] = plus_var
        entry["orf_plus_widgets"] = add_file_row(card, 6, plus_label, plus_tip, plus_var,
                                                 button_text="Browse  ORF + flanks",
                                                 button_width=32)

        source_var.trace_add("write", lambda *_, e=entry: apply_orf_flank_source(e))

        sel_var = tk.StringVar(value=PLACEHOLDER)
        entry["codon_selection_file_path"] = sel_var
        add_file_row(card, 7, "Codon selection (optional)", "OPTIONAL. Choose specific chimera "
                     "insertion points for THIS ORF from an .xlsx file (residue numbers in "
                     "column 1); overrides the codon sampling frequency below. Leave unset to "
                     "scan by the sampling frequency.", sel_var,
                     button_text="Codon Selection: by File (optional)", button_width=36,
                     optional=True)

        # Graphical alternative to the .xlsx: pick insertion codons from the protein.
        entry["codon_positions"] = []
        picked_var = tk.StringVar(value=NO_CODONS_PICKED)
        entry["_picked_var"] = picked_var
        picked_lbl = ttk.Label(card, textvariable=picked_var, style="Path.TLabel",
                               wraplength=560, anchor="w", justify="left")
        picked_lbl.grid(row=8, column=0, sticky="we", padx=(14, 6), pady=10)
        pick_btn = ttk.Button(card, text="Codon Selection: by Codon Picker (optional)",
                              style="Browse.TButton", width=36,
                              command=lambda e=entry: open_codon_picker(e))
        pick_btn.grid(row=8, column=1, sticky="e", padx=(6, 14), pady=10)
        attach_tip(picked_lbl, "OPTIONAL. Pick specific insertion codons graphically from this ORF's "
                   "protein sequence (translated from the ORF FASTA). An alternative to the "
                   ".xlsx above; if both are set, the picked codons are added to it.", pick_btn)

        orf_entries.append(entry)
        renumber_orfs()
        apply_flank_mode_to_entry(entry)
        return entry

    # Folder-button help text; the per-ORF version names the flank files, the
    # global version makes clear those flanks are no longer required per ORF.
    folder_tip_per_orf = (
        "Discover ORFs from a folder. Files are named '<gene>_coding.fa', "
        "'<gene>_flank5.fa', '<gene>_flank3.fa' (+ optional '<gene>_codonSelection.xlsx'). "
        "The gene name is the part before the suffix.")
    folder_tip_global = (
        "Discover ORFs from a folder. Files are named '<gene>_coding.fa' (+ optional "
        "'<gene>_codonSelection.xlsx'); the gene name is the part before the suffix. The "
        "global 5'/3' flanks above are applied to every ORF, so per-ORF flank files are "
        "not needed.")
    folder_tip = None   # assigned once the folder button is created, below

    # Remembered on the first refresh, while the holder is still packed: the
    # sibling it sits above, so re-showing it restores the original order
    # instead of dropping it at the bottom of the form.
    flank_holder_anchor = []

    def refresh_flank_mode(*_):
        glob = flank_mode.get() == "global"
        if not flank_holder_anchor:
            slaves = content.pack_slaves()
            i = slaves.index(global_flank_holder)
            flank_holder_anchor.append(slaves[i + 1] if i + 1 < len(slaves) else None)
        if glob:
            anchor = flank_holder_anchor[0]
            if anchor is not None and anchor.winfo_exists():
                global_flank_holder.pack(fill="x", before=anchor)
            else:
                global_flank_holder.pack(fill="x")
        else:
            global_flank_holder.pack_forget()
        for entry in orf_entries:
            apply_flank_mode_to_entry(entry)
        if folder_tip is not None:
            folder_tip.text = folder_tip_global if glob else folder_tip_per_orf

    flank_mode.trace_add("write", refresh_flank_mode)

    # Buttons: add one ORF, or discover many from a folder.
    orf_buttons = ttk.Frame(content, style="TFrame")
    orf_buttons.pack(anchor="w", pady=(0, 4))
    ttk.Button(orf_buttons, text="+  Add ORF", style="Browse.TButton", command=add_orf).pack(
        side="left", padx=(0, 8))

    def populate_orfs(orfs):
        """Replace the queued ORF cards with the discovered ORFs."""
        for entry in list(orf_entries):
            entry["card"].destroy()
        orf_entries.clear()
        for orf in orfs:
            entry = add_orf()
            entry["geneName"].set(orf.get("geneName", ""))
            for key in ("orf_file_path", "flank5_file_path", "flank3_file_path",
                        "codon_selection_file_path"):
                if orf.get(key):
                    entry[key].set(orf[key])
        renumber_orfs()

    def protein_fastas_in(folder):
        """Return (name, path, accession) for each FASTA in *folder* that is protein."""
        found = []
        for name in sorted(os.listdir(folder)):
            if os.path.splitext(name)[1].lower() not in _FASTA_EXTS:
                continue
            path = os.path.join(folder, name)
            try:
                if fetch_cds.fasta_is_protein(path):
                    found.append((name, path, fetch_cds.accession_from_fasta(path)))
            except OSError:
                continue
        return found

    def fetch_cds_then_load(folder, proteins):
        """Fetch a CDS from UniProt for each protein file, then re-discover the folder."""
        status_var.set("Fetching %d coding sequence(s) from UniProt…" % len(proteins))
        folder_btn.config(state="disabled")
        clear_console()
        console_banner("===== Fetching %d coding sequence(s) from UniProt =====\n" % len(proteins))

        def worker():
            written, failed = [], []
            for i, (name, _path, accession) in enumerate(proteins):
                if not accession:
                    console_error("  %s: no UniProt accession in header\n" % name)
                    failed.append("%s (no UniProt accession in header)" % name)
                    continue
                console_log("  %s (%s): fetching…\n" % (name, accession))
                try:
                    result = fetch_cds.fetch_cds_for_accession(accession)
                    fetch_cds.write_coding_fasta(folder, result)
                    console_log("    -> %s  %s [%d bp]\n"
                                % (result.gene, result.refseq_id, len(result.sequence)))
                    written.append(result.gene)
                except fetch_cds.ServiceUnavailable as exc:
                    # Service is down: stop rather than retry every remaining gene.
                    remaining = len(proteins) - i - 1
                    console_error("    %s\n" % exc)
                    console_error("  Aborting: UniProt/NCBI is unavailable; %d gene(s) not "
                                  "attempted.\n" % remaining)
                    failed.append("%s (service unavailable)" % name)
                    break
                except Exception as exc:   # per-file: report and keep going
                    console_error("    failed: %s\n" % exc)
                    failed.append("%s (%s)" % (name, exc))
                if i < len(proteins) - 1:
                    time.sleep(0.34)       # respect NCBI's rate limit
            root.after(0, lambda: finish_fetch(folder, written, failed))

        threading.Thread(target=worker, daemon=True).start()

    def finish_fetch(folder, written, failed):
        folder_btn.config(state="normal")
        try:
            orfs, _skipped = discover_orf_folder(folder)
        except OSError as exc:
            messagebox.showerror("Folder error", str(exc))
            return
        populate_orfs(orfs)
        msg = "Fetched %d coding sequence(s); loaded %d ORF(s)." % (len(written), len(orfs))
        if failed:
            msg += "\nCould not fetch: %s" % "; ".join(failed)
        status_var.set(msg)
        messagebox.showinfo("UniProt fetch complete", msg)

    def load_folder():
        folder = filedialog.askdirectory(initialdir=dialog_dir["path"])
        if not folder:
            return
        remember_dir(folder)
        try:
            orfs, skipped = discover_orf_folder(folder)
        except OSError as exc:
            messagebox.showerror("Folder error", str(exc))
            return

        # Any protein FASTAs present can't be scanned; offer to fetch their CDS.
        proteins = protein_fastas_in(folder)
        if proteins and not orfs:
            names = ", ".join(name for name, _, _ in proteins)
            if messagebox.askyesno(
                    "Protein sequences found",
                    "%d file(s) here are protein sequences, not DNA:\n%s\n\n"
                    "PAM-scanning needs the coding DNA (ATG..stop). Fetch each gene's CDS "
                    "from UniProt (via its RefSeq cross-reference) and use those instead?\n\n"
                    "Requires an internet connection." % (len(proteins), names)):
                fetch_cds_then_load(folder, proteins)
                return

        if not orfs:
            extra = ""
            if proteins:
                extra = ("\n\n(%d protein FASTA(s) were found but not fetched.)" % len(proteins))
            messagebox.showwarning(
                "No ORFs found",
                "No ORF (DNA) files were found in:\n%s%s" % (folder, extra))
            return

        populate_orfs(orfs)
        msg = "Loaded %d ORF(s) from folder." % len(orfs)
        if skipped:
            msg += "  Ignored %d unrecognized file(s): %s" % (len(skipped), ", ".join(skipped))
        status_var.set(msg)
        if skipped:
            messagebox.showinfo("Folder loaded", msg)

    folder_btn = ttk.Button(orf_buttons, text="Load folder…", style="Browse.TButton",
                            command=load_folder)
    folder_btn.pack(side="left")
    folder_tip = attach_tip(folder_btn, folder_tip_per_orf)

    # --- Shared input files ---------------------------------------------
    card = section("Shared input files", 0)
    for r, (key, blabel, tip) in enumerate(SHARED_FILE_FIELDS):
        # A field whose own label says "(optional)" says so on its path line too.
        add_file_row(card, r, blabel, tip, shared_file_vars[key],
                     optional="(optional)" in blabel)

    # Local BLAST database: browse to any member file; store the -db prefix.
    db_row = len(SHARED_FILE_FIELDS)
    db_lbl = ttk.Label(card, textvariable=blast_db_var, style="Path.TLabel",
                       wraplength=560, anchor="w", justify="left")
    db_lbl.grid(row=db_row, column=0, sticky="we", padx=(14, 6), pady=10)

    def browse_blast_db():
        chosen = filedialog.askopenfilename(
            initialdir=dialog_dir["path"],
            title="Select any file of an existing BLAST database (e.g. yeast.nin)")
        if chosen:
            remember_dir(chosen)
            blast_db_var.set(blast_db_prefix(chosen))

    db_btns = ttk.Frame(card, style="Card.TFrame")
    db_btns.grid(row=db_row, column=1, sticky="e", padx=(6, 14), pady=10)
    ttk.Button(db_btns, text="Use genome (auto)", style="Small.TButton",
               command=lambda: blast_db_var.set(AUTO_DB)).pack(side="left", padx=(0, 6))
    db_btn = ttk.Button(db_btns, text="Browse database…", style="Browse.TButton",
                        command=browse_blast_db, width=20)
    db_btn.pack(side="left")
    attach_tip(db_lbl,
               "The BLAST+ database for off-target checks. By default it is built once from the "
               "genome above (cached in ~/.pam_scanning/blastdb) and reused, so you don't need to "
               "build one yourself. Optionally Browse to an existing prebuilt database (any member "
               "file, e.g. yeast.nin) to use that instead; 'Use genome (auto)' returns to the default.",
               db_btn)

    # --- Sequence & primer settings -------------------------------------
    card = section("Sequence & primer settings", 1)
    for r, (key, label, _default, tip) in enumerate(STRING_FIELDS):
        add_entry_row(card, r, label, tip, string_vars[key], mono=True)

    # --- Scan parameters -------------------------------------------------
    card = section("Scan parameters", 1)
    for r, (key, label, _default, tip) in enumerate(INT_FIELDS):
        add_entry_row(card, r, label, tip, int_vars[key], validate=True)

    # --- Output ----------------------------------------------------------
    card = section("Output", 0)
    out_lbl = ttk.Label(card, textvariable=output_var, style="Path.TLabel",
                        wraplength=560, anchor="w", justify="left")
    out_lbl.grid(row=0, column=0, sticky="we", padx=(14, 6), pady=10)

    def browse_dir():
        chosen = filedialog.askdirectory(initialdir=dialog_dir["path"])
        if chosen:
            remember_dir(chosen)
            output_var.set(chosen)

    out_btn = ttk.Button(card, text="Choose directory", style="Browse.TButton",
                         command=browse_dir, width=26)
    out_btn.grid(row=0, column=1, sticky="e", padx=(6, 14), pady=10)
    attach_tip(out_lbl, "Directory where the time-stamped PAM-scan output folder(s) are written.", out_btn)

    # --- Action bar ------------------------------------------------------
    action = ttk.Frame(content, style="TFrame")
    action.pack(fill="x", pady=(22, 6))
    run_button = ttk.Button(action, text="▶  Run PAM Scan", style="Accent.TButton")
    run_button.pack()
    ttk.Label(content, textvariable=status_var, style="Status.TLabel").pack(anchor="center", pady=(8, 4))

    # Start with one ORF and the initial (per-ORF) flank mode applied.
    add_orf()
    refresh_flank_mode()

    # --- Run logic -------------------------------------------------------
    def collect_shared():
        shared = {key: var.get() for key, var in shared_file_vars.items()}
        shared.update({key: var.get() for key, var in string_vars.items()})
        db_text = blast_db_var.get()
        # AUTO_DB (the default) => empty, so pamscan builds/caches a db from the genome.
        shared["localBlastDb"] = "" if db_text == AUTO_DB else blast_db_prefix(db_text)
        if flank_mode.get() == "global":
            # Send exactly one kwarg per flank, so pamscan never sees file+sequence.
            for key, _label, short, _tip in GLOBAL_FLANK_FIELDS:
                if global_flank_source[key].get() == "file":
                    shared[key + "_file_path"] = global_flank_path_vars[key].get()
                    continue
                text = global_flank_seq_vars[key].get().strip()
                if not text:
                    shared[key + "_sequence"] = ""   # validate() reports the empty box
                    continue
                try:
                    shared[key + "_sequence"] = parse_sequence_text(text, "%s sequence" % short)
                except ValueError as exc:
                    messagebox.showerror("Invalid sequence", str(exc))
                    return None
        for key, var in int_vars.items():
            text = var.get().strip()
            if not text:
                messagebox.showerror("Invalid input", "Please enter a value for '%s'." % key)
                return None
            shared[key] = int(text)
        out = output_var.get().strip()
        shared["outputPath"] = out if out and out != PLACEHOLDER else "."
        return shared

    def collect_orfs():
        per_orf_mode = flank_mode.get() == "per_orf"
        orfs = []
        for entry in orf_entries:
            orf_path = entry["orf_file_path"].get()
            gene = entry["geneName"].get().strip()
            # Blank gene name: derive it from the ORF file name.
            if not gene and _is_set(orf_path) and os.path.isfile(orf_path):
                gene = gene_name_from_orf_path(orf_path)
            orf = {
                "geneName": gene,
                "orf_file_path": orf_path,
                "codon_selection_file_path": entry["codon_selection_file_path"].get(),
            }
            positions = entry.get("codon_positions") or []
            if positions:
                orf["codon_selection_positions"] = list(positions)
            if per_orf_mode:
                if entry["flank_source"].get() == "combined":
                    orf["orf_plus_file_path"] = entry["orf_plus_file_path"].get()
                else:
                    orf["flank5_file_path"] = entry["flank5_file_path"].get()
                    orf["flank3_file_path"] = entry["flank3_file_path"].get()
            orfs.append(orf)
        return orfs

    def _is_set(value):
        return bool(value) and value != PLACEHOLDER

    def validate(shared, orfs):
        if not _is_set(shared["local_genome_file_path"]):
            messagebox.showerror("Missing input", "Please select the genome sequence.")
            return False
        if not orfs:
            messagebox.showerror("Missing input", "Please add at least one ORF.")
            return False
        if flank_mode.get() == "global":
            for key, _label, short, _tip in GLOBAL_FLANK_FIELDS:
                if global_flank_source[key].get() == "file":
                    if not _is_set(shared.get(key + "_file_path")):
                        messagebox.showerror(
                            "Missing input",
                            "Global flank mode: please select the %s FASTA file." % short)
                        return False
                elif not shared.get(key + "_sequence"):
                    messagebox.showerror(
                        "Missing input",
                        "Global flank mode: please enter the %s sequence." % short)
                    return False
        for i, orf in enumerate(orfs, start=1):
            if not orf["geneName"].strip():
                messagebox.showerror("Missing input",
                                     "ORF %d: enter a gene name, or select an ORF file to "
                                     "derive one from." % i)
                return False
            if not _is_set(orf["orf_file_path"]):
                messagebox.showerror("Missing input",
                                     "ORF %d (%s): please select the ORF file." % (i, orf["geneName"]))
                return False
            if flank_mode.get() == "per_orf":
                if "orf_plus_file_path" in orf:
                    needed = (("orf_plus_file_path", "ORF + flanks"),)
                else:
                    needed = (("flank5_file_path", "5' flank"),
                              ("flank3_file_path", "3' flank"))
                for key, name in needed:
                    if not _is_set(orf[key]):
                        messagebox.showerror("Missing input",
                                             "ORF %d (%s): please select the %s file."
                                             % (i, orf["geneName"], name))
                        return False
        return True

    def finish(error, out_path, n):
        run_button.config(state="normal")
        if error is not None:
            status_var.set("Error — see the progress console.")
            console_error("\nError: %s\n" % error)   # no popup; surfaced in the console
        else:
            status_var.set("Done. %d ORF(s) written under: %s" % (n, out_path))

    def start_scan_worker(shared, orfs):
        out_path = shared["outputPath"]
        n = len(orfs)

        def worker():
            writer = _ConsoleWriter(console_log)
            try:
                from pam_scanning.chimeras import pamscan
                with contextlib.redirect_stdout(writer):
                    for i, orf in enumerate(orfs, start=1):
                        console_banner("\n===== ORF %d/%d: %s =====\n" % (i, n, orf["geneName"]))
                        msg = "Running ORF %d/%d (%s)… see the progress console." % (
                            i, n, orf["geneName"])
                        root.after(0, lambda m=msg: status_var.set(m))
                        result = pamscan(**dict(shared, **orf))
                        if isinstance(result, dict) and result.get("plot_png"):
                            root.after(0, lambda p=result["plot_png"]: embed_plot(p))
                console_banner("\n===== Done: %d ORF(s) written =====\n" % n)
                root.after(0, lambda: finish(None, out_path, n))
            except Exception as exc:  # surface any failure back on the UI thread
                console_error("\nError: %s\n" % exc)
                root.after(0, lambda e=exc: finish(e, out_path, n))

        threading.Thread(target=worker, daemon=True).start()

    def blast_install_failed(exc):
        run_button.config(state="normal")
        status_var.set("BLAST+ install failed — see console.")
        console_error("\n%s\n" % exc)
        messagebox.showerror("BLAST+ install failed", str(exc))

    def ensure_blast_then_scan(shared, orfs):
        """Run the scan; if BLAST+ is missing, offer to install it first."""
        from pam_scanning import blast_setup

        if blast_setup.ensure_available() is not None:
            start_scan_worker(shared, orfs)
            return
        if not messagebox.askyesno(
                "BLAST+ not found",
                "BLAST+ (blastn) is required for off-target screening but was not found on "
                "this computer.\n\nDownload it now automatically? The official NCBI BLAST+ "
                "binaries are placed in ~/.pam_scanning/blast — nothing is added to conda or "
                "your environment. This may take a few minutes; progress is shown in the "
                "console."):
            run_button.config(state="normal")
            status_var.set("BLAST+ is required to run a scan.")
            return
        status_var.set("Downloading BLAST+…")
        console_banner("===== Installing BLAST+ (NCBI) =====\n")

        def worker():
            try:
                blast_setup.install_blast(log=console_log)
                root.after(0, lambda: start_scan_worker(shared, orfs))
            except Exception as exc:
                root.after(0, lambda e=exc: blast_install_failed(e))

        threading.Thread(target=worker, daemon=True).start()

    def run_scan():
        shared = collect_shared()
        if shared is None:
            return
        orfs = collect_orfs()
        if not validate(shared, orfs):
            return
        run_button.config(state="disabled")
        clear_console()
        ensure_blast_then_scan(shared, orfs)

    run_button.config(command=run_scan)

    # Bias the initial split toward the form; the sash stays user-draggable. The
    # wider window gives the progress console (and its embedded plot) more room.
    root.after(60, lambda: paned.sashpos(0, 860))

    root.mainloop()


if __name__ == "__main__":
    main()
