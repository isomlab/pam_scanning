# PAM-scanning

Design CRISPR/Cas9 guide RNAs and chimera-insertion primers **across an ORF** — every
codon by default, or a chosen subset ("PAM scanning"). Given an open reading frame and
its genomic context, the tool:

1. finds all NGG/CCN PAM sites across the ORF and the guides that target them,
2. **silences** each PAM site (or, failing that, the guide body) with synonymous
   codon substitutions so re-cutting is prevented after editing,
3. uses **NCBI BLAST+** to screen guides for off-target potential against the host
   genome (including conservative "PAM-inclusion" tracking),
4. selects the optimal guide for each insertion codon, and
5. emits ready-to-order guide and insertion primers (96-/384-well plate layout),
   QC FASTA files, a scannable-sequence report, and off-target warnings.

This is the reference implementation accompanying the PAM-scanning manuscript.

For a step-by-step account of the algorithm — how PAM sites become silenced,
off-target-screened, per-codon-optimal guides and orderable primers — see
[`docs/pipeline.md`](docs/pipeline.md).

## Contents

- [Installation](#installation) — conda, pip, or a double-click launcher
- [Prerequisites: a local BLAST database](#prerequisites-a-local-blast-database)
- [Usage](#usage) — [command line](#command-line) and [graphical interface](#graphical-interface)
- [Inputs](#inputs)
- **Docs:** [pipeline overview](docs/pipeline.md) · [full usage reference](docs/usage.md) · [BLAST setup](docs/blast_setup.md) · [non-coder install guide](INSTALL.md)

## Installation

PAM-scanning depends on the external **NCBI BLAST+** toolkit (`blastn`), which is most
easily obtained through conda. The recommended path is therefore conda.

### Easiest (for lab members / non-coders)

Install [Miniforge](https://conda-forge.org/download/) once, then double-click the
launcher for your computer in [`launchers/`](launchers):
**`PAM Scanning.command`** (Mac) or **`PAM Scanning.bat`** (Windows). The first run
sets everything up automatically (Python, the app, and BLAST+); after that it just
opens the graphical app. Step-by-step guide: [`INSTALL.md`](INSTALL.md).

### conda (recommended)

```bash
git clone https://github.com/isomlab/pam_scanning.git
cd pam_scanning
conda env create -f environment.yml
conda activate pam_scanning
```

This installs Python, `openpyxl`, **BLAST+**, and the `pam_scanning` package itself
(editable). Once a [Bioconda](https://bioconda.github.io/) release is published you will
also be able to run:

```bash
conda install -c bioconda -c conda-forge pam_scanning
```

### pip

```bash
pip install pam_scanning
```

`pip` does **not** install BLAST+ — install it separately (e.g.
`conda install -c bioconda blast`, or from the
[NCBI BLAST+ executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)).
If `blastn` is missing when you run a scan, the GUI offers to download BLAST+ for you
(and the CLI accepts `--install-blast`) — the official NCBI binaries are placed in
`~/.pam_scanning/blast`, with no conda required and nothing added to your environment.

## Prerequisites: a local BLAST database

You need a local BLAST+ nucleotide database for your host genome. See
[`docs/blast_setup.md`](docs/blast_setup.md) for step-by-step instructions
(`makeblastdb`, where the genome FASTA goes, and how the database name maps to the
`--blast-db` option).

## Usage

### Command line

```bash
pam-scan \
    --orf examples/fasta/S288C_YBL016W_FUS3_coding.fa \
    --flank5 examples/fasta/S288C_YBL016W_FUS3_flank5.fa \
    --flank3 examples/fasta/S288C_YBL016W_FUS3_flank3.fa \
    --gene-name Fus3 \
    --output ./results
```

No `--genome` or `--blast-db` is needed for a standard yeast scan: the bundled *S. cerevisiae*
BY4741 genome is used by default, and its BLAST database is built once (with `makeblastdb`) and
cached automatically. To scan against a different yeast species/strain/variant, pass your own
`--genome`; give `--blast-db` only to reuse a prebuilt database.

To scan **multiple ORFs** in one run, either list them in a tab-separated manifest
(`--manifest examples/manifest.tsv`) or point at a folder of conventionally-named FASTA
files (`--orf-dir examples/orf_folder`), sharing the same genome.
Flanks can be given per ORF or shared globally via a single `--flank5`/`--flank3` pair;
see [`docs/usage.md`](docs/usage.md#multiple-orfs).

Working from **UniProt**? Those downloads are protein sequences, not coding DNA. The
bundled `pam-scan-fetch-cds` command fetches each entry's CDS (via its RefSeq
cross-reference, choosing the isoform that matches the canonical protein) and writes
`‹gene›_coding.fa` ready for `--orf-dir`; see
[`docs/usage.md`](docs/usage.md#preparing-orfs-from-uniprot-pam-scan-fetch-cds).

All parameters can also be supplied via `--config run.toml` (or `run.json`); explicit
flags override config values. Run `pam-scan --help` for the full list. The yeast codon
table is bundled and used by default; supply `--codon-table` to override it. To restrict
insertion to specific residues, give either `--codon-selection` (an `.xlsx`) or
`--codon-positions "52, 89, 100-105"`. See [`docs/usage.md`](docs/usage.md) for the full
parameter reference and a description of the output directory.

### Graphical interface

```bash
pam-scan-gui
```

A Tkinter form that collects the same parameters and runs the identical pipeline.
Use **+ Add ORF** to queue ORFs one at a time, or **Load folder…** to discover a folder
of them; the *Flank inputs* control switches between per-ORF and global 5′/3′ flanks.
In global mode each flank can be loaded **From file** or typed in with **Enter sequence**.
Each ORF offers two ways to restrict insertion to specific residues:
**Codon Selection: by File** (an `.xlsx`) and **Codon Selection: by Codon Picker**.
The picker opens the gene's translated protein as a numbered, clickable grid — click
a residue to add or remove it, **drag to select a run**, or type positions/ranges.
Every pick appears in a list below, grouped into runs, with a **Remove** button on
each line so a selection can be pruned before running.

## Inputs

| Input | Description |
| --- | --- |
| ORF FASTA (`--orf`) | The coding sequence, ATG → stop. |
| 5′ flank (`--flank5` / `--flank5-seq`) | The 100 bp immediately upstream of the ATG (the `-` side), so positions at the start of the ORF can be scanned. Give a FASTA file or a literal sequence. |
| 3′ flank (`--flank3` / `--flank3-seq`) | The 100 bp immediately downstream of the stop (the `+` side), so positions at the end of the ORF can be scanned. Give a FASTA file or a literal sequence. |
| Genome FASTA (`--genome`) | The **yeast** host genome for off-target evaluation (also the source for your BLAST DB). PAM scanning is always run in yeast — the ORF is ported in from its source organism — so this is always a yeast genome. Defaults to the bundled BY4741 genome; override to pick a different yeast species/strain/variant. |
| Codon table | Codon-usage table; the bundled yeast table is used if omitted. |
| Codon selection (optional) | Specific residues to target, as an `.xlsx` (`--codon-selection`), a position list (`--codon-positions`), or the GUI **Pick codons…** picker. Omit to scan every codon (or every *N*th via the sampling gap). |

For a batch of ORFs, a TSV manifest (`--manifest`) supplies the ORF, 5′/3′ flank, and
optional codon-selection paths one row at a time; see [`docs/usage.md`](docs/usage.md#multiple-orfs).

## Development

```bash
pip install -e ".[test]"
pytest
```

The unit tests run without BLAST+. The end-to-end regression test
(`tests/test_end_to_end.py`) is skipped automatically unless `blastn`/`makeblastdb`
are on your `PATH`.

## Citation

If you use this software, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff).

## License

[MIT](LICENSE) © 2026 Daniel Isom (Isom Lab)
