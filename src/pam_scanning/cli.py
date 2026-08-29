"""Command-line interface for PAM-scanning.

This is a thin, scriptable front-end over :func:`pam_scanning.chimeras.pamscan`.
It exposes exactly the parameters the GUI collects, so GUI and CLI runs are
equivalent and reproducible. Parameters may be given as flags and/or loaded from
a ``--config`` file (JSON or TOML); explicit flags override config values.

Examples
--------
Scan a single ORF against a local yeast BLAST database::

    pam-scan \\
        --orf examples/fasta/S288C_YBL016W_FUS3_coding.fa \\
        --flank5 examples/fasta/S288C_YBL016W_FUS3_flank5.fa \\
        --flank3 examples/fasta/S288C_YBL016W_FUS3_flank3.fa \\
        --genome /path/to/BY4741_Toronto_2012.fsa \\
        --blast-db yeast \\
        --gene-name Fus3 \\
        --output ./results

Scan many ORFs from a tab-separated manifest (one ORF per row), with the shared
parameters supplied by flags or ``--config``::

    pam-scan --manifest orfs.tsv --genome genome.fsa --blast-db yeast --output ./results

The manifest has a header row and one row per ORF. Recognized columns:
``gene``, ``orf``, ``flank5``, ``flank3`` (required) and ``codon_selection``
(optional). Relative paths are resolved against the manifest's directory.

Either flank may instead be given as a literal sequence with ``--flank5-seq`` /
``--flank3-seq``, which is convenient for the global flanks shared by a whole
batch of ORFs::

    pam-scan --orf-dir ./orfs --flank5-seq "$(cat up.txt)" --flank3-seq ACGT... \\
        --genome genome.fsa --blast-db yeast

For a given side, pass the FASTA flag or the sequence flag, not both. A per-ORF
flank from a manifest or folder takes precedence over a global flank sequence.

Instead of two separate flank files, the ORF and both its flanks may be given as
ONE FASTA with ``--orf-plus``::

    pam-scan --orf FUS3_coding.fa --orf-plus FUS3_orfPlus.fa \\
        --gene-name Fus3 --genome genome.fsa --output ./results

``--orf`` is still required: the ORF is located inside the combined sequence by
exact match and the flanks are taken from either side, so the flanks may be any
length rather than a fixed 100 bp. The run fails with a clear message if the ORF
is absent from that file, occurs more than once, or leaves a flank of zero
length. ``--orf-plus`` is mutually exclusive with ``--flank5``/``--flank3`` and
their ``-seq`` forms. In a manifest the column is ``orf_plus``; in a folder the
filename suffix is ``_orfPlus`` (also ``_withFlanks``).
"""

import argparse
import csv
import os
import re
import shutil
import sys

from pam_scanning.chimeras import parse_sequence_text, parse_codon_positions


# A BLAST database is a set of files sharing a prefix; 'blastn -db' takes that
# prefix. This matches a trailing BLAST member-file extension (e.g. '.nin',
# '.psq', or a multi-volume '.00.nhr') so a selected member file can be mapped
# back to the prefix.
_BLAST_DB_SUFFIX = re.compile(r"\.(?:\d{2,4}\.)?[pn][a-z]{2}$", re.IGNORECASE)


def blast_db_prefix(value):
    """Return the ``blastn -db`` prefix for a database name or member-file path.

    If *value* ends in a BLAST database extension (``.nin``, ``.psq``, a
    multi-volume ``.00.nhr``, etc.) the extension is stripped so the result is
    the path/name passed to ``blastn -db``. A bare name like ``yeast`` (resolved
    via ``$BLASTDB``) is returned unchanged.
    """
    if not value:
        return value
    return _BLAST_DB_SUFFIX.sub("", value)


# Parameter defaults shared with the GUI. Keys match the kwargs that
# pam_scanning.chimeras.pamscan expects.
DEFAULTS = {
    "geneName": None,
    "localBlastDb": "",   # empty => auto-build a database from the genome (cached)
    "guidePrimerForwardSuffix": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAG",
    "insertPrimerForwardSuffix": "GAAGATGTTGTCTGTTGCTCTATGTCATAT",
    "insertPrimerReverseSuffix": "CTTCTACAACAGACAACGAGATACAGTATA",
    "primerLength": 100,
    "maxPamCutGap": 60,
    "codonsSamplingGap": 1,
    "pamInclusionThreshold": 5,
    "pamInclusionSequenceThreshold": 15,
}

# Keys that vary per ORF (supplied by single-ORF flags or by manifest rows).
PER_ORF_KEYS = ("geneName", "orf_file_path", "flank5_file_path", "flank3_file_path",
                "orf_plus_file_path", "codon_selection_file_path")

# chimeras.pamscan treats this sentinel as "not provided".
NOT_SELECTED = "No file selected"


def _is_set(value):
    """True when an option holds a real value rather than being unset/sentinel."""
    return bool(value) and value != NOT_SELECTED


def _build_parser():
    p = argparse.ArgumentParser(
        prog="pam-scan",
        description="Design PAM-scanning chimera-insertion guides and primers for one or more ORFs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Inputs / outputs.
    p.add_argument("--config", help="JSON or TOML file with any of the options below.")
    p.add_argument("--manifest",
                   help="TSV of ORFs (columns: gene, orf, flank5, flank3, [codon_selection]); "
                        "one ORF per row. Shared params come from flags/--config.")
    p.add_argument("--orf-dir", dest="orf_dir",
                   help="Folder of ORFs named '<gene>_coding.fa', '<gene>_flank5.fa', "
                        "'<gene>_flank3.fa' (+ optional '<gene>_codonSelection.xlsx'). "
                        "Per-ORF flanks in the folder are used; for GLOBAL flanks shared by "
                        "every ORF, omit them and pass --flank5/--flank3 instead.")
    p.add_argument("--orf", dest="orf_file_path", help="ORF FASTA (ATG..stop). [single ORF]")
    p.add_argument("--flank5", dest="flank5_file_path",
                   help="5' flank FASTA: 100 bp upstream of the ATG (the '-' side). Single ORF, "
                        "or a GLOBAL 5' flank applied to every ORF in --manifest/--orf-dir.")
    p.add_argument("--flank3", dest="flank3_file_path",
                   help="3' flank FASTA: 100 bp downstream of the stop (the '+' side). Single ORF, "
                        "or a GLOBAL 3' flank applied to every ORF in --manifest/--orf-dir.")
    p.add_argument("--orf-plus", dest="orf_plus_file_path",
                   help="FASTA holding the 5' flank, the ORF and the 3' flank as ONE "
                        "sequence. An alternative to --flank5/--flank3: the ORF is located "
                        "inside it by exact match and the flanks are taken from either side, "
                        "so any flank length works. --orf is still required. Mutually "
                        "exclusive with --flank5/--flank3 and their -seq forms.")
    p.add_argument("--flank5-seq", dest="flank5_sequence",
                   help="5' flank as a literal sequence instead of a FASTA file. Mutually "
                        "exclusive with --flank5.")
    p.add_argument("--flank3-seq", dest="flank3_sequence",
                   help="3' flank as a literal sequence instead of a FASTA file. Mutually "
                        "exclusive with --flank3.")
    p.add_argument("--genome", dest="local_genome_file_path",
                   help="Yeast host genome FASTA for off-target checks. PAM scanning is always "
                        "run in yeast (the ORF is ported in), so this is always a yeast genome; "
                        "use it to choose the yeast species/strain/variant. Defaults to the "
                        "bundled BY4741 genome when omitted.")
    p.add_argument("--codon-table", dest="codon_table_file_path",
                   help="Codon-usage table. Defaults to the bundled yeast table.")
    p.add_argument("--codon-selection", dest="codon_selection_file_path",
                   help="Optional .xlsx of specific insertion sites (overrides sampling). [single ORF]")
    p.add_argument("--codon-positions", dest="codon_selection_positions",
                   type=parse_codon_positions,
                   help="Specific insertion codons as positions/ranges, e.g. '52, 89, 100-105' "
                        "(overrides sampling; adds to --codon-selection if both are given).")
    p.add_argument("-o", "--output", dest="outputPath", default=".",
                   help="Directory to write the time-stamped results into.")
    p.add_argument("--install-blast", dest="install_blast", action="store_true",
                   help="If BLAST+ ('blastn') is missing, download the official NCBI binaries "
                        "into ~/.pam_scanning/blast (no conda needed) instead of exiting.")

    # String parameters.
    p.add_argument("--gene-name", dest="geneName", help="Label used in output filenames. [single ORF]")
    p.add_argument("--blast-db", dest="localBlastDb",
                   help="Optional prebuilt BLAST+ database (name or path). When omitted, a "
                        "database is built once from --genome and cached in ~/.pam_scanning/blastdb.")
    p.add_argument("--guide-primer-forward-suffix", dest="guidePrimerForwardSuffix")
    p.add_argument("--insert-primer-forward-suffix", dest="insertPrimerForwardSuffix")
    p.add_argument("--insert-primer-reverse-suffix", dest="insertPrimerReverseSuffix")

    # Integer parameters.
    p.add_argument("--primer-length", dest="primerLength", type=int)
    p.add_argument("--max-pam-cut-gap", dest="maxPamCutGap", type=int)
    p.add_argument("--codon-sampling-gap", dest="codonsSamplingGap", type=int)
    p.add_argument("--max-pam-inclusions", dest="pamInclusionThreshold", type=int)
    p.add_argument("--max-pam-inclusion-length", dest="pamInclusionSequenceThreshold", type=int)

    return p


def _load_config(path):
    """Load a JSON or TOML config file into a dict of pamscan kwargs."""
    ext = os.path.splitext(path)[1].lower()
    with open(path, "rb") as fh:
        if ext == ".json":
            import json

            return json.load(fh)
        if ext in (".toml", ".tml"):
            try:
                import tomllib  # Python >= 3.11
            except ModuleNotFoundError:  # pragma: no cover
                try:
                    import tomli as tomllib
                except ModuleNotFoundError:
                    raise SystemExit(
                        "Reading TOML config requires Python 3.11+ or the 'tomli' package."
                    )
            return tomllib.load(fh)
    raise SystemExit("Config file must end in .json or .toml: %s" % path)


# Manifest column name -> pamscan kwarg key. Several spellings are accepted.
_MANIFEST_COLUMNS = {
    "gene": "geneName", "genename": "geneName", "name": "geneName",
    "orf": "orf_file_path",
    "flank5": "flank5_file_path", "5flank": "flank5_file_path", "upstream": "flank5_file_path",
    "flank3": "flank3_file_path", "3flank": "flank3_file_path", "downstream": "flank3_file_path",
    "orf_plus": "orf_plus_file_path", "orfplus": "orf_plus_file_path",
    "orf+flanks": "orf_plus_file_path", "orf_flanks": "orf_plus_file_path",
    "codon_selection": "codon_selection_file_path",
    "codon-selection": "codon_selection_file_path",
    "codonselection": "codon_selection_file_path",
}


def _load_manifest(path):
    """Parse a TSV manifest into a list of per-ORF kwarg dicts.

    Relative file paths in the manifest are resolved against the manifest's own
    directory so a manifest plus its FASTA files can live together.
    """
    base = os.path.dirname(os.path.abspath(path))

    def resolve(value):
        value = (value or "").strip()
        if not value:
            return None
        return value if os.path.isabs(value) else os.path.normpath(os.path.join(base, value))

    file_keys = {"orf_file_path", "flank5_file_path", "flank3_file_path",
                 "codon_selection_file_path"}
    orfs = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("Manifest is empty: %s" % path)
        for raw in reader:
            entry = {}
            for col, value in raw.items():
                key = _MANIFEST_COLUMNS.get((col or "").strip().lower())
                if key is None:
                    continue
                value = (value or "").strip()
                if not value:
                    continue
                entry[key] = resolve(value) if key in file_keys else value
            if not entry:
                continue  # blank line
            orfs.append(entry)
    if not orfs:
        raise SystemExit("Manifest has no ORF rows: %s" % path)
    return orfs


# --- Folder discovery (flat folder, filename convention) -------------------

# Recognized extensions for the two kinds of input file.
_FASTA_EXTS = (".fa", ".fasta", ".fna", ".fas")
_XLSX_EXTS = (".xlsx",)

# Filename-suffix -> role (pamscan kwarg). The gene name is the filename stem
# with the matched suffix removed; matching is case-insensitive. Flank and
# codon-selection suffixes are checked before the ORF suffixes (no stem ends in
# more than one, but this keeps intent explicit).
_ROLE_SUFFIXES = (
    ("orf_plus_file_path",
     ("_orfplus", "_orf_plus", "_orfcontext", "_withflanks", "_plusflanks")),
    ("flank5_file_path", ("_flank5", "_5flank", "_upstream", "_5prime", "_5p")),
    ("flank3_file_path", ("_flank3", "_3flank", "_downstream", "_3prime", "_3p")),
    ("codon_selection_file_path",
     ("_codonselection", "_codon_selection", "_codonselect", "_codons")),
    ("orf_file_path", ("_coding", "_orf", "_cds")),
)

# A trailing RefSeq nucleotide/protein accession in a derived gene name (e.g.
# 'ABCB1_NM_001348945.2' from 'ABCB1_NM_001348945.2_ORF.fasta'). Stripped so the
# gene name is just the symbol. Names without such a token are left unchanged, so
# the '<strain>_<systematic>_<symbol>' convention is unaffected.
_REFSEQ_ACCESSION = re.compile(r"_(?:NM|XM|NR|XR|NP|XP)_\d+(?:\.\d+)?$", re.IGNORECASE)


def _gene_symbol(name):
    """Trim a trailing RefSeq accession from a derived gene name."""
    return _REFSEQ_ACCESSION.sub("", name)


# ORF filename suffixes, longest first, so '_coding'/'_orf'/'_cds' are stripped.
_ORF_NAME_SUFFIXES = ("_coding", "_orf", "_cds")


def gene_name_from_orf_path(path):
    """Derive a gene name from an ORF FASTA filename (same rule as folder discovery).

    Strips the extension, a trailing ORF role suffix (``_coding`` / ``_orf`` / ``_cds``),
    and any trailing RefSeq accession, so ``AKT1_NM_001382431.1_ORF.fasta`` yields
    ``AKT1`` and ``FUS3_coding.fa`` yields ``FUS3``.
    """
    stem = os.path.splitext(os.path.basename(path))[0]
    low = stem.lower()
    for suffix in _ORF_NAME_SUFFIXES:
        if low.endswith(suffix):
            stem = stem[: len(stem) - len(suffix)]
            break
    return _gene_symbol(stem)


def discover_orf_folder(path):
    """Discover ORFs in a flat folder by filename convention (deterministic).

    Files are named ``<gene><suffix><ext>``; the suffix marks the role:

    * ORF:           ``_coding`` / ``_orf`` / ``_cds``          (FASTA)
    * 5' flank:      ``_flank5`` / ``_5flank`` / ``_upstream``  (FASTA)
    * 3' flank:      ``_flank3`` / ``_3flank`` / ``_downstream``(FASTA)
    * ORF + flanks:  ``_orfPlus`` / ``_withFlanks``              (FASTA)
    * codon select.: ``_codonSelection`` / ``_codons``          (.xlsx)

    The gene name is the filename stem with the role suffix removed, and a
    trailing RefSeq accession trimmed (so 'ABCB1_NM_001348945.2_ORF.fasta' yields
    'ABCB1'). Files are grouped into one ORF per gene. Per-ORF flanks found in the
    folder are kept; when absent, the caller can supply a global 5'/3' flank instead.

    Returns ``(orfs, skipped)`` where *orfs* is a list of per-ORF kwarg dicts
    ordered by gene name, and *skipped* is a list of FASTA/.xlsx file names whose
    suffix matched no known role (so callers can warn rather than silently drop).
    """
    grouped = {}   # geneName -> {role: full path}
    skipped = []
    for name in sorted(os.listdir(path)):
        full = os.path.join(path, name)
        if not os.path.isfile(full):
            continue
        stem, ext = os.path.splitext(name)
        ext = ext.lower()
        if ext not in _FASTA_EXTS and ext not in _XLSX_EXTS:
            continue  # not an input file we care about (READMEs, .tsv, etc.)
        low = stem.lower()
        role = gene = None
        for key, suffixes in _ROLE_SUFFIXES:
            match = next((s for s in suffixes if low.endswith(s)), None)
            if match is not None:
                role = key
                gene = _gene_symbol(stem[: len(stem) - len(match)])
                break
        if role is None:
            skipped.append(name)
            continue
        # The extension must suit the role (codon selection is a spreadsheet).
        wants_xlsx = role == "codon_selection_file_path"
        if wants_xlsx != (ext in _XLSX_EXTS):
            skipped.append(name)
            continue
        grouped.setdefault(gene, {"geneName": gene})[role] = full

    orfs = []
    for gene in sorted(grouped):
        entry = grouped[gene]
        if "orf_file_path" in entry:
            orfs.append(entry)
        else:
            # A gene with flanks but no ORF file is not a scannable ORF.
            skipped.extend(os.path.basename(v) for k, v in entry.items() if k != "geneName")
    return orfs, skipped


def _check_blast(install=False):
    """Ensure the external NCBI BLAST+ 'blastn' executable is available.

    With *install* set, BLAST+ is installed via conda/bioconda if it is missing;
    otherwise a missing 'blastn' exits with instructions (including how to let the
    tool install it).
    """
    from pam_scanning import blast_setup

    if blast_setup.ensure_available() is not None:
        return
    if install:
        try:
            blast_setup.install_blast(log=lambda text: print(text, end=""))
        except RuntimeError as exc:
            sys.exit("Error: %s" % exc)
        return
    sys.exit(
        "Error: 'blastn' (NCBI BLAST+) was not found on your PATH.\n"
        "PAM-scanning requires BLAST+ for off-target evaluation.\n"
        "Let this tool download it for you:     pam-scan ... --install-blast\n"
        "or install it yourself with conda:     conda install -c bioconda blast\n"
        "or see https://www.ncbi.nlm.nih.gov/books/NBK279690/"
    )


def build_kwargs(argv=None):
    """Resolve config file + command-line flags into a (base kwargs, args) pair.

    For a single ORF the base is the complete kwargs dict; for a ``--manifest``
    run it is the shared base that each ORF row is merged onto.
    """
    args = _build_parser().parse_args(argv)

    kwargs = dict(DEFAULTS)
    if args.config:
        kwargs.update(_load_config(args.config))

    # Flags explicitly provided on the command line override config/defaults.
    for key, value in vars(args).items():
        if key in ("config", "manifest", "orf_dir", "install_blast") or value is None:
            continue
        kwargs[key] = value

    # File-path defaults: chimeras.pamscan treats these sentinels as "not provided".
    kwargs.setdefault("codon_table_file_path", NOT_SELECTED)
    kwargs.setdefault("codon_selection_file_path", NOT_SELECTED)
    kwargs.setdefault("outputPath", ".")
    # Accept either a database name or a path to one of its member files.
    kwargs["localBlastDb"] = blast_db_prefix(kwargs.get("localBlastDb"))
    _normalize_flank_sequences(kwargs)
    return kwargs, args


def _normalize_flank_sequences(kwargs):
    """Validate and clean any flank given as a sequence; reject file+sequence clashes.

    Parsing here means a typo in a pasted sequence fails immediately, rather than
    after the run has already reached BLAST.
    """
    for side in ("5", "3"):
        seq_key, path_key = "flank%s_sequence" % side, "flank%s_file_path" % side
        raw = kwargs.get(seq_key)
        if not _is_set(raw):
            kwargs.pop(seq_key, None)
            continue
        if _is_set(kwargs.get(path_key)):
            sys.exit("Error: use either --flank%s or --flank%s-seq, not both." % (side, side))
        try:
            kwargs[seq_key] = parse_sequence_text(raw, "%s' flank sequence" % side)
        except ValueError as exc:
            sys.exit("Error: %s" % exc)


def _validate(kwargs):
    """Fail early, with clear messages, on missing required inputs for one ORF."""
    required = {
        "orf_file_path": "--orf / manifest 'orf'",
        "local_genome_file_path": "--genome",
    }
    missing = [flag for key, flag in required.items() if not _is_set(kwargs.get(key))]
    if _is_set(kwargs.get("orf_plus_file_path")):
        # One file carries both flanks, so the separate flank inputs must be absent.
        clashes = [flag for key, flag in (
            ("flank5_file_path", "--flank5 / manifest 'flank5'"),
            ("flank3_file_path", "--flank3 / manifest 'flank3'"),
            ("flank5_sequence", "--flank5-seq"),
            ("flank3_sequence", "--flank3-seq"),
        ) if _is_set(kwargs.get(key))]
        if clashes:
            sys.exit("Error: --orf-plus supplies both flanks, so do not also give: %s"
                     % ", ".join(clashes))
    else:
        # Each flank may arrive as a FASTA file or as a sequence; one of the two is required.
        for side in ("5", "3"):
            if not _is_set(kwargs.get("flank%s_file_path" % side)) \
                    and not _is_set(kwargs.get("flank%s_sequence" % side)):
                missing.append("--flank%s / --flank%s-seq / --orf-plus / manifest 'flank%s'"
                               % (side, side, side))
    if missing:
        sys.exit("Error: missing required input(s): %s" % ", ".join(sorted(missing)))

    if not kwargs.get("geneName"):
        sys.exit("Error: a gene name is required (--gene-name or manifest 'gene').")

    # Only the file form has a path to check; sequences were validated at parse time.
    for key, flag in (("orf_file_path", "--orf"),
                      ("orf_plus_file_path", "--orf-plus"),
                      ("flank5_file_path", "--flank5"),
                      ("flank3_file_path", "--flank3"),
                      ("local_genome_file_path", "--genome")):
        path = kwargs.get(key)
        if not _is_set(path):
            continue
        if not os.path.isfile(path):
            sys.exit("Error: file for %s does not exist: %s" % (flag, path))


def _merge_orf(base, orf):
    """Merge one ORF's fields onto the shared base.

    A flank supplied for this ORF (manifest column or folder file) wins over the
    global flank, whether that global flank was a file or an entered sequence --
    so the two forms can never both reach pamscan for the same side.
    """
    kwargs = dict(base)
    kwargs.update(orf)
    for side in ("5", "3"):
        if _is_set(orf.get("flank%s_file_path" % side)):
            kwargs.pop("flank%s_sequence" % side, None)
    return kwargs


def main(argv=None):
    """Console-script entry point for ``pam-scan``.

    Handles a single ORF, a ``--manifest`` TSV, or an ``--orf-dir`` folder. In the
    batch modes, any 5'/3' flank not given per ORF falls back to a global
    ``--flank5``/``--flank3`` flag (or ``--flank5-seq``/``--flank3-seq``), so
    global and per-ORF flanks both work.
    """
    base, args = build_kwargs(argv)
    if args.manifest and args.orf_dir:
        sys.exit("Error: use either --manifest or --orf-dir, not both.")
    _check_blast(install=args.install_blast)

    # Default the off-target genome to the bundled yeast genome when none is given.
    if not _is_set(base.get("local_genome_file_path")):
        from pam_scanning.chimeras import default_genome_path
        base["local_genome_file_path"] = default_genome_path()

    from pam_scanning.chimeras import pamscan

    orfs = None
    if args.manifest:
        orfs = _load_manifest(args.manifest)
    elif args.orf_dir:
        orfs, skipped = discover_orf_folder(args.orf_dir)
        if skipped:
            print("Warning: ignored %d unrecognized file(s) in %s: %s"
                  % (len(skipped), args.orf_dir, ", ".join(skipped)), file=sys.stderr)
        if not orfs:
            sys.exit("Error: no ORF files (e.g. '<gene>_coding.fa') found in: %s" % args.orf_dir)

    if orfs is not None:
        total = len(orfs)
        for i, orf in enumerate(orfs, start=1):
            kwargs = _merge_orf(base, orf)           # per-ORF fields override the shared base
            _validate(kwargs)
            print("=== PAM scan %d/%d: %s ===" % (i, total, kwargs.get("geneName")))
            pamscan(**kwargs)
        return 0

    _validate(base)
    pamscan(**base)
    return 0


if __name__ == "__main__":
    sys.exit(main())
