# Usage reference

PAM-scanning can be driven from the command line (`pam-scan`), a Tkinter GUI
(`pam-scan-gui`), or programmatically via `pam_scanning.chimeras.pamscan(**kwargs)`.
All three share the same parameters and produce the same output.

For *what* the tool computes (the algorithm behind these parameters), see
[`pipeline.md`](pipeline.md); this page is the parameter and output reference.

## Parameters

| CLI flag | kwarg | Default | Description |
| --- | --- | --- | --- |
| `--orf` | `orf_file_path` | *(required, per ORF)* | ORF FASTA, ATG → stop. |
| `--flank5` | `flank5_file_path` | *(one flank per side required)* | FASTA of the 100 bp immediately **upstream** of the ATG (the `-` side). Lets the scan reach positions at the start of the ORF. |
| `--flank5-seq` | `flank5_sequence` | *(one flank per side required)* | The 5′ flank as a literal sequence instead of a FASTA file. Mutually exclusive with `--flank5`. [†](#notes-on-specific-options) |
| `--flank3` | `flank3_file_path` | *(one flank per side required)* | FASTA of the 100 bp immediately **downstream** of the stop (the `+` side). Lets the scan reach positions at the end of the ORF. |
| `--flank3-seq` | `flank3_sequence` | *(one flank per side required)* | The 3′ flank as a literal sequence instead of a FASTA file. Mutually exclusive with `--flank3`. |
| `--orf-plus` | `orf_plus_file_path` | *(none)* | FASTA holding the 5′ flank, the ORF and the 3′ flank as **one** sequence, instead of `--flank5`/`--flank3`. `--orf` is still required and locates the ORF inside it, so the flanks may be any length. Mutually exclusive with `--flank5`/`--flank3` and their `-seq` forms. |
| `--manifest` | *(n/a)* | *(none)* | TSV of ORFs (one row each) for batch runs; see [Multiple ORFs](#multiple-orfs). |
| `--genome` | `local_genome_file_path` | bundled BY4741 genome | The **yeast** host genome FASTA for off-target checks. [†](#notes-on-specific-options) |
| `--blast-db` | `localBlastDb` | *(auto)* | **Optional.** A prebuilt BLAST+ database. Built automatically from `--genome` when omitted. [†](#notes-on-specific-options) |
| `--gene-name` | `geneName` | *(required, per ORF)* | Label used in output filenames. |
| `--codon-table` | `codon_table_file_path` | bundled yeast table | Codon-usage table (`.cusp`-style). |
| `--codon-selection` | `codon_selection_file_path` | *(none, per ORF)* | `.xlsx` of specific residues to target; overrides sampling. |
| `--codon-positions` | `codon_selection_positions` | *(none)* | Specific insertion codons as a position/range list, e.g. `"52, 89, 100-105"` (1-based). [†](#notes-on-specific-options) |
| `--output` / `-o` | `outputPath` | `.` | Directory to write the time-stamped run into. |
| `--guide-primer-forward-suffix` | `guidePrimerForwardSuffix` | `GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAG` | Suffix that amplifies the CRISPR plasmid after the targeting sequence. |
| `--insert-primer-forward-suffix` | `insertPrimerForwardSuffix` | `GAAGATGTTGTCTGTTGCTCTATGTCATAT` | 5′→3′ insertion-primer suffix (chimera payload, forward). |
| `--insert-primer-reverse-suffix` | `insertPrimerReverseSuffix` | `CTTCTACAACAGACAACGAGATACAGTATA` | 3′→5′ insertion-primer suffix (chimera payload, reverse). |
| `--primer-length` | `primerLength` | `100` | Total length (bp) of the chimera-amplification primers. |
| `--max-pam-cut-gap` | `maxPamCutGap` | `60` | Max bp between sequential PAM cut sites (the empirical 30 bp rule × 2). |
| `--codon-sampling-gap` | `codonsSamplingGap` | `1` | Insert at every Nth codon (1 = exhaustive). Ignored when a codon-selection file is given. |
| `--max-pam-inclusions` | `pamInclusionThreshold` | `5` | Max allowed PAM inclusions per guide solution. |
| `--max-pam-inclusion-length` | `pamInclusionSequenceThreshold` | `15` | Min matched length (bp) counted as a PAM inclusion. |

### Notes on specific options

Marked **†** in the table above.

- **`--flank5-seq` / `--flank3-seq`** — accept A/C/G/T/N; a pasted FASTA record or base
  numbering is tolerated and stripped.

- **`--genome`** — PAM scanning is always performed in yeast (the ORF is ported in from
  its source organism), so this is always a yeast genome. The bundled *S. cerevisiae*
  BY4741 genome ships gzipped and is expanded to `~/.pam_scanning/genome` on first use.
  Override it to use a different yeast species, strain, or variant.

- **`--blast-db`** — when omitted, a BLAST+ database is built once from `--genome` with
  `makeblastdb`, cached in `~/.pam_scanning/blastdb` (keyed to the genome), and reused,
  so you never have to build one by hand. To override, give a bare name (resolved via
  `$BLASTDB`), a path prefix, or a path to any member file (e.g. `/data/yeast.nin`,
  which is reduced to the prefix). In the GUI, **Browse database…** overrides and
  **Use genome (auto)** restores the default.

- **`--codon-positions`** — overrides sampling, and adds to `--codon-selection` if both
  are given. In the GUI this is the **Pick codons…** picker.

### Config file

Instead of flags you can pass `--config run.toml` (or `.json`). Keys are the **kwarg**
names from the table above. Flags given on the command line override config values.

```toml
# run.toml
orf_file_path = "examples/fasta/S288C_YBL016W_FUS3_coding.fa"
flank5_file_path = "examples/fasta/S288C_YBL016W_FUS3_flank5.fa"
flank3_file_path = "examples/fasta/S288C_YBL016W_FUS3_flank3.fa"
local_genome_file_path = "/path/to/BY4741_Toronto_2012.fsa"
localBlastDb = "yeast"
geneName = "Fus3"
outputPath = "./results"
codonsSamplingGap = 1
```

```bash
pam-scan --config run.toml
```

## Choosing which codons to target

By default the chimera is inserted at **every** codon (or every *N*th, via
`--codon-sampling-gap`). To restrict insertion to specific residues instead, use
any one of:

- **`--codon-positions "52, 89, 100-105"`** — a 1-based list of positions and
  inclusive ranges. Codon *n* is residue *n* of the protein (codon 1 = the start Met).
- **`--codon-selection sites.xlsx`** — a spreadsheet with a `Mutations` sheet whose
  first column holds the residue numbers (columns 2–3 may name the original/target
  residue but are not required for site selection).
- **GUI → Codon Selection: by Codon Picker** — each ORF card has a button that opens
  the gene's translated protein as a numbered, clickable grid. Click a residue to add
  or remove it, **drag to select a run**, or type the same `"52, 89, 100-105"` syntax.
  Picked codons appear in a list below the grid — contiguous runs are grouped onto one
  line — each with its own **Remove** button, so a selection can be pruned before
  running. The picked codons are shown on the ORF card and passed to the run; it is
  the graphical equivalent of `--codon-positions`. (The **Codon Selection: by File**
  button on the same card is the `.xlsx` route.)

Any of these overrides the sampling gap. If both `--codon-positions` (or the picker)
and `--codon-selection` are supplied, their residues are combined.

When a specific set is requested, the **scannability grid** (`<gene>-scannableMap.pdf/png`)
reports only those codons: unselected codons are greyed out, while a selected codon that
is still inaccessible stays white so it stands out.

## Multiple ORFs

Several ORFs can be scanned in one invocation; each produces its own time-stamped
run directory. There are two ways to supply them.

### Manifest (TSV)

List one ORF per row in a tab-separated **manifest** and supply the shared
parameters with flags or `--config`. Recognized columns: `gene`, `orf`, `flank5`,
`flank3` (or `orf_plus` in their place) and `codon_selection` (optional); relative paths resolve
against the manifest's directory. See `examples/manifest.tsv`.

```tsv
gene	orf	flank5	flank3	codon_selection
Fus3	fasta/FUS3_coding.fa	fasta/FUS3_flank5.fa	fasta/FUS3_flank3.fa	codon_selection/Fus3.xlsx
Kss1	fasta/KSS1_coding.fa	fasta/KSS1_flank5.fa	fasta/KSS1_flank3.fa
```

```bash
pam-scan --manifest examples/manifest.tsv \
    --genome /path/to/BY4741_Toronto_2012.fsa --blast-db yeast --output ./results
```

### Folder (`--orf-dir`)

Point at a folder of FASTA files named by convention; the gene name is the part
**before** the role suffix. See `examples/orf_folder/`.

| Role | Suffix (any of) | Type |
| --- | --- | --- |
| ORF | `_coding`, `_orf`, `_cds` | FASTA |
| 5′ flank | `_flank5`, `_5flank`, `_upstream` | FASTA |
| 3′ flank | `_flank3`, `_3flank`, `_downstream` | FASTA |
| codon selection (optional) | `_codonSelection`, `_codons` | `.xlsx` |

```bash
pam-scan --orf-dir examples/orf_folder \
    --genome /path/to/BY4741_Toronto_2012.fsa --blast-db yeast --output ./results
```

Files whose suffix matches no role are ignored with a warning.

### Global vs per-ORF flanks

Flanks can be supplied **per ORF** (each ORF row/file has its own 5′/3′ flank) or
**globally** (one 5′/3′ pair applied to every ORF). For global flanks, omit the
per-ORF flank columns/files and pass a global flank once; it fills in for any ORF
that doesn't specify its own. A global flank can be either a FASTA file
(`--flank5`/`--flank3`) or a literal sequence (`--flank5-seq`/`--flank3-seq`):

```bash
# Global flanks (from files) shared by every ORF in the folder:
pam-scan --orf-dir examples/orf_folder \
    --flank5 shared_flank5.fa --flank3 shared_flank3.fa \
    --genome genome.fsa --blast-db yeast --output ./results

# The same, giving the flanks as sequences instead of files:
pam-scan --orf-dir examples/orf_folder \
    --flank5-seq "$(cat shared_flank5.txt)" --flank3-seq ACGT...100bp \
    --genome genome.fsa --blast-db yeast --output ./results
```

A per-ORF flank (from a manifest row or a folder file) always takes precedence
over a global flank for that ORF and side.

In the GUI, use **+ Add ORF** to queue ORFs one at a time or **Load folder…** to
discover them, and choose **Per-ORF flanks** or **Global flanks** under *Flank inputs*.
In global mode each flank has a **From file** / **Enter sequence** toggle; the
sequence box shows a live base count (and reads *invalid* on a non-DNA character).

## Preparing ORFs from UniProt (`pam-scan-fetch-cds`)

PAM-scanning needs the **coding DNA** (ATG→stop), but UniProt distributes
**protein** FASTA. The `pam-scan-fetch-cds` helper bridges the two so a folder of
UniProt downloads becomes a folder of PAM-scannable ORFs. The source organism only
supplies the ORF — the scan itself is always performed in yeast, so `--genome`
(off-target evaluation) remains a yeast genome regardless of where the gene came from.

```bash
# From a folder of UniProt protein FASTAs (accession read from each header):
pam-scan-fetch-cds --protein-dir ./uniprot_fastas --email you@example.com

# Or from explicit accessions, written to a chosen folder:
pam-scan-fetch-cds P60709 P08183 P31749 -o ./orfs --email you@example.com
```

Each protein is resolved deterministically:

1. the UniProt entry's **RefSeq** cross-reference gives the curated mRNA(s);
2. when an entry has several mRNA isoforms, the one whose CDS **translates to the
   UniProt canonical protein** is chosen (not an arbitrary first), so the exact
   ORF is never guessed at;
3. NCBI E-utilities returns that mRNA's CDS, written as `‹gene›_coding.fa` — the
   name [`--orf-dir`](#folder---orf-dir) discovery recognizes.

The stop codon is retained. Each file's header records provenance
(`CDS from RefSeq NM_… (UniProt …)`), and any anomaly (no start ATG, no stop,
frame, unmatched isoform) is reported per accession rather than silently written.
`--email` is passed to NCBI as a courtesy and is recommended for large batches;
`--delay` (default 0.34 s) keeps within NCBI's rate limit.

In the **GUI**, this is automatic: if **Load folder…** finds protein sequences
instead of DNA, it offers to fetch their CDS from UniProt (the same routine) and
loads the results — so a folder of UniProt downloads can be used directly.

> The helper writes `‹gene›_coding.fa` alongside (or into `-o`); the original
> protein FASTAs, having no role suffix, are simply ignored by folder discovery.
> Use **Global flanks** (a shared 5′/3′ pair) since these files carry no per-ORF
> flanks. An internet connection is required (UniProt + NCBI); it is the only part
> of PAM-scanning that goes online.

## Output

Each run creates a time-stamped directory `‹gene›-chimera-insertions-‹YYYY.MM.DD-HH.MM.SS›/`
under `--output`, containing:

| Path | Contents |
| --- | --- |
| `QC/‹gene›-guideSolutions.xlsx` | Per-codon optimal guide, silenced guide, cut gap, PAM inclusions, and insertion primers. |
| `QC/‹gene›-scannableSequence.txt` | Fraction of the ORF that is PAM-scannable, plus the masked sequence and affected codons. |
| `QC/solutionsFasta/‹gene›-‹codon›.fa` | SnapGene-viewable silenced ORF for each insertion site. |
| `QC/` (copies) | The input ORF and 5′/3′ flank FASTAs, plus the assembled `‹gene›-orfPlusContext.fa`, kept for provenance. |
| `ORDER/‹gene›-primerOrder.xlsx` | Plate-laid-out guide + insertion primer order (96- or 384-well). |
| `BLAST+/` | Raw and final `blastn` query/result files for off-target review. |
| `WARNINGS/‹gene›-pamInclusionWarnings*.txt` | Guides carrying potential PAM inclusions (conservative + super-conservative). |
| `WARNINGS/‹gene›-unscannableCodons.txt` | Codons for which no acceptable guide could be found. |

A flank supplied as a literal sequence (`--flank5-seq` / `--flank3-seq`) is written out
as `‹gene›-flank5.fa` / `‹gene›-flank3.fa`, so a run is reproducible from its own
`QC/` folder either way.
