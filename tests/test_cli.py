"""Unit tests for the command-line front-end (no BLAST required).

These cover the new flank-based inputs and the multi-ORF manifest, exercising
the argument/config/manifest plumbing without invoking the BLAST pipeline.
"""

import os

import pytest

from pam_scanning import cli


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTA = os.path.join(REPO_ROOT, "examples", "fasta")
ORF = os.path.join(FASTA, "S288C_YBL016W_FUS3_coding.fa")
FLANK5 = os.path.join(FASTA, "S288C_YBL016W_FUS3_flank5.fa")
FLANK3 = os.path.join(FASTA, "S288C_YBL016W_FUS3_flank3.fa")


def test_single_orf_flags_build_kwargs():
    base, args = cli.build_kwargs([
        "--orf", ORF, "--flank5", FLANK5, "--flank3", FLANK3,
        "--genome", "/tmp/genome.fsa", "--gene-name", "Fus3",
    ])
    assert args.manifest is None
    assert base["orf_file_path"] == ORF
    assert base["flank5_file_path"] == FLANK5
    assert base["flank3_file_path"] == FLANK3
    assert base["geneName"] == "Fus3"
    # Untouched defaults survive.
    assert base["primerLength"] == 100
    assert base["localBlastDb"] == ""   # empty => auto-build a db from the genome


def test_validate_requires_both_flanks():
    base, _ = cli.build_kwargs([
        "--orf", ORF, "--flank5", FLANK5,
        "--genome", ORF, "--gene-name", "Fus3",
    ])
    with pytest.raises(SystemExit):
        cli._validate(base)  # missing --flank3


def test_validate_passes_for_complete_single_orf():
    base, _ = cli.build_kwargs([
        "--orf", ORF, "--flank5", FLANK5, "--flank3", FLANK3,
        "--genome", ORF, "--gene-name", "Fus3",
    ])
    cli._validate(base)  # should not raise


def test_flank_sequence_flag_is_parsed_and_normalized():
    base, _ = cli.build_kwargs([
        "--orf", ORF, "--flank5-seq", "> hdr\nac gt\nNNN", "--flank3", FLANK3,
        "--genome", ORF, "--gene-name", "Fus3",
    ])
    # Header, whitespace stripped; upper-cased; file-path form left unset.
    assert base["flank5_sequence"] == "ACGTNNN"
    assert base.get("flank5_file_path") is None
    cli._validate(base)  # a sequence flank satisfies the requirement


def test_flank_file_and_sequence_are_mutually_exclusive():
    with pytest.raises(SystemExit):
        cli.build_kwargs([
            "--orf", ORF, "--flank5", FLANK5, "--flank5-seq", "ACGT",
            "--flank3", FLANK3, "--genome", ORF, "--gene-name", "Fus3",
        ])


def test_invalid_flank_sequence_exits_early():
    with pytest.raises(SystemExit):
        cli.build_kwargs([
            "--orf", ORF, "--flank5-seq", "ACGTX", "--flank3", FLANK3,
            "--genome", ORF, "--gene-name", "Fus3",
        ])


def test_validate_requires_a_flank_in_each_direction():
    base, _ = cli.build_kwargs([
        "--orf", ORF, "--flank5-seq", "ACGT",
        "--genome", ORF, "--gene-name", "Fus3",
    ])
    with pytest.raises(SystemExit):
        cli._validate(base)  # no 3' flank, file or sequence


def test_per_orf_flank_file_overrides_global_sequence():
    base = {"flank5_sequence": "AAAA", "flank3_sequence": "TTTT", "geneName": "shared"}
    merged = cli._merge_orf(base, {"flank5_file_path": "/data/GeneA_flank5.fa"})
    assert merged["flank5_file_path"] == "/data/GeneA_flank5.fa"
    assert "flank5_sequence" not in merged      # file wins; no clash reaches pamscan
    assert merged["flank3_sequence"] == "TTTT"  # untouched side keeps the global sequence


def test_manifest_parses_rows_and_resolves_paths(tmp_path):
    manifest = tmp_path / "orfs.tsv"
    # Relative paths must resolve against the manifest's own directory.
    (tmp_path / "a.fa").write_text(">a\nATG\n")
    (tmp_path / "u.fa").write_text(">u\nAAA\n")
    (tmp_path / "d.fa").write_text(">d\nTTT\n")
    manifest.write_text(
        "gene\torf\tflank5\tflank3\n"
        "GeneA\ta.fa\tu.fa\td.fa\n"
        "GeneB\ta.fa\tu.fa\td.fa\n"
    )
    orfs = cli._load_manifest(str(manifest))
    assert [o["geneName"] for o in orfs] == ["GeneA", "GeneB"]
    assert orfs[0]["orf_file_path"] == str(tmp_path / "a.fa")
    assert orfs[0]["flank5_file_path"] == str(tmp_path / "u.fa")
    assert orfs[0]["flank3_file_path"] == str(tmp_path / "d.fa")


def test_manifest_accepts_column_aliases(tmp_path):
    manifest = tmp_path / "orfs.tsv"
    manifest.write_text(
        "name\torf\tupstream\tdownstream\tcodon_selection\n"
        "G1\to.fa\tup.fa\tdn.fa\tsel.xlsx\n"
    )
    orfs = cli._load_manifest(str(manifest))
    assert orfs[0]["geneName"] == "G1"
    assert orfs[0]["flank5_file_path"].endswith("up.fa")
    assert orfs[0]["flank3_file_path"].endswith("dn.fa")
    assert orfs[0]["codon_selection_file_path"].endswith("sel.xlsx")


def test_bundled_manifest_example_loads():
    example = os.path.join(REPO_ROOT, "examples", "manifest.tsv")
    orfs = cli._load_manifest(example)
    assert len(orfs) == 1
    assert orfs[0]["geneName"] == "Fus3"
    assert os.path.isfile(orfs[0]["orf_file_path"])
    assert os.path.isfile(orfs[0]["flank5_file_path"])
    assert os.path.isfile(orfs[0]["flank3_file_path"])


# --- Folder discovery (--orf-dir) ------------------------------------------

def _touch(path, text=">x\nACGT\n"):
    with open(path, "w") as fh:
        fh.write(text)


def test_discover_folder_groups_by_gene_and_role(tmp_path):
    _touch(tmp_path / "FUS3_coding.fa")
    _touch(tmp_path / "FUS3_flank5.fa")
    _touch(tmp_path / "FUS3_flank3.fa")
    (tmp_path / "FUS3_codonSelection.xlsx").write_bytes(b"PK\x03\x04stub")
    _touch(tmp_path / "KSS1_orf.fasta")
    _touch(tmp_path / "KSS1_upstream.fa")
    _touch(tmp_path / "KSS1_downstream.fa")

    orfs, skipped = cli.discover_orf_folder(str(tmp_path))
    assert skipped == []
    assert [o["geneName"] for o in orfs] == ["FUS3", "KSS1"]
    fus3 = orfs[0]
    assert fus3["orf_file_path"].endswith("FUS3_coding.fa")
    assert fus3["flank5_file_path"].endswith("FUS3_flank5.fa")
    assert fus3["flank3_file_path"].endswith("FUS3_flank3.fa")
    assert fus3["codon_selection_file_path"].endswith("FUS3_codonSelection.xlsx")
    # Alias suffixes resolve to the same roles.
    kss1 = orfs[1]
    assert kss1["orf_file_path"].endswith("KSS1_orf.fasta")
    assert kss1["flank5_file_path"].endswith("KSS1_upstream.fa")
    assert kss1["flank3_file_path"].endswith("KSS1_downstream.fa")


def test_discover_folder_trims_refseq_accession_to_symbol(tmp_path):
    _touch(tmp_path / "ABCB1_NM_001348945.2_ORF.fasta")
    _touch(tmp_path / "ACTB_NM_001101.5_ORF.fasta")
    orfs, skipped = cli.discover_orf_folder(str(tmp_path))
    assert skipped == []
    assert [o["geneName"] for o in orfs] == ["ABCB1", "ACTB"]


def test_discover_folder_keeps_symbol_last_names(tmp_path):
    # No RefSeq token: the '<strain>_<systematic>_<symbol>' name is left intact.
    _touch(tmp_path / "S288C_YBL016W_FUS3_coding.fa")
    orfs, _skipped = cli.discover_orf_folder(str(tmp_path))
    assert [o["geneName"] for o in orfs] == ["S288C_YBL016W_FUS3"]


def test_gene_symbol_strips_only_refseq_accessions():
    assert cli._gene_symbol("ABCB1_NM_001348945.2") == "ABCB1"
    assert cli._gene_symbol("AKT1_XM_017001.1") == "AKT1"
    assert cli._gene_symbol("FUS3") == "FUS3"                 # nothing to strip
    assert cli._gene_symbol("S288C_YBL016W_FUS3") == "S288C_YBL016W_FUS3"


def test_discover_folder_global_flank_layout(tmp_path):
    """Only ORF files present: flanks are expected to come from a global flag."""
    _touch(tmp_path / "GeneA_coding.fa")
    _touch(tmp_path / "GeneB_coding.fa")
    orfs, skipped = cli.discover_orf_folder(str(tmp_path))
    assert skipped == []
    assert [o["geneName"] for o in orfs] == ["GeneA", "GeneB"]
    assert all("flank5_file_path" not in o for o in orfs)


def test_discover_folder_reports_unrecognized_and_orphans(tmp_path):
    _touch(tmp_path / "GeneA_coding.fa")
    _touch(tmp_path / "mystery.fa")          # no role suffix
    _touch(tmp_path / "GeneC_flank5.fa")     # flank with no ORF -> orphan
    orfs, skipped = cli.discover_orf_folder(str(tmp_path))
    assert [o["geneName"] for o in orfs] == ["GeneA"]
    assert set(skipped) == {"mystery.fa", "GeneC_flank5.fa"}


def test_bundled_orf_folder_example_loads():
    folder = os.path.join(REPO_ROOT, "examples", "orf_folder")
    orfs, skipped = cli.discover_orf_folder(folder)
    assert skipped == []
    assert len(orfs) == 1
    o = orfs[0]
    assert o["geneName"] == "FUS3"
    assert os.path.isfile(o["orf_file_path"])
    assert os.path.isfile(o["flank5_file_path"])
    assert os.path.isfile(o["flank3_file_path"])
    assert os.path.isfile(o["codon_selection_file_path"])


def test_gene_name_from_orf_path():
    assert cli.gene_name_from_orf_path("/x/FUS3_coding.fa") == "FUS3"
    assert cli.gene_name_from_orf_path("/x/KSS1_orf.fna") == "KSS1"
    assert cli.gene_name_from_orf_path("/x/AKT1_NM_001382431.1_ORF.fasta") == "AKT1"
    # No recognized suffix / a strain-systematic-symbol name is left intact.
    assert cli.gene_name_from_orf_path("/x/S288C_YBL016W_FUS3_coding.fa") == "S288C_YBL016W_FUS3"
    assert cli.gene_name_from_orf_path("/x/plainname.fasta") == "plainname"


def test_manifest_and_orf_dir_are_mutually_exclusive():
    with pytest.raises(SystemExit):
        cli.main(["--manifest", "m.tsv", "--orf-dir", "d", "--genome", "g.fsa"])


# --- BLAST availability check ----------------------------------------------

def test_check_blast_exits_when_missing_without_install(monkeypatch):
    from pam_scanning import blast_setup
    monkeypatch.setattr(blast_setup, "ensure_available", lambda: None)
    with pytest.raises(SystemExit, match="install-blast"):
        cli._check_blast(install=False)


def test_check_blast_passes_when_present(monkeypatch):
    from pam_scanning import blast_setup
    monkeypatch.setattr(blast_setup, "ensure_available", lambda: "/usr/bin/blastn")
    cli._check_blast(install=False)   # should not raise


def test_check_blast_installs_when_requested(monkeypatch):
    from pam_scanning import blast_setup
    calls = {"installed": False}
    monkeypatch.setattr(blast_setup, "ensure_available", lambda: None)

    def fake_install(log):
        calls["installed"] = True
        return "/home/u/.pam_scanning/blast/ncbi-blast/bin/blastn"

    monkeypatch.setattr(blast_setup, "install_blast", fake_install)
    cli._check_blast(install=True)
    assert calls["installed"] is True


# --- BLAST database prefix -------------------------------------------------

def test_blast_db_prefix_strips_member_extensions():
    assert cli.blast_db_prefix("/data/yeast.nin") == "/data/yeast"
    assert cli.blast_db_prefix("/data/yeast.nsq") == "/data/yeast"
    assert cli.blast_db_prefix("/data/prot.psq") == "/data/prot"
    # Multi-volume member file.
    assert cli.blast_db_prefix("/data/big.00.nhr") == "/data/big"


def test_blast_db_prefix_leaves_names_and_unrelated_paths_alone():
    assert cli.blast_db_prefix("yeast") == "yeast"
    assert cli.blast_db_prefix("/data/yeast") == "/data/yeast"
    assert cli.blast_db_prefix("/genomes/genome.fsa") == "/genomes/genome.fsa"
    assert cli.blast_db_prefix("") == ""


def test_build_kwargs_normalizes_blast_db_path():
    base, _ = cli.build_kwargs(["--blast-db", "/data/yeast.nin", "--orf", ORF])
    assert base["localBlastDb"] == "/data/yeast"


def test_codon_positions_flag_parses_into_kwarg():
    base, _ = cli.build_kwargs([
        "--orf", ORF, "--flank5", FLANK5, "--flank3", FLANK3,
        "--genome", "/tmp/genome.fsa", "--gene-name", "Fus3",
        "--codon-positions", "10, 20-22",
    ])
    assert base["codon_selection_positions"] == [10, 20, 21, 22]


def test_codon_positions_absent_by_default():
    base, _ = cli.build_kwargs([
        "--orf", ORF, "--flank5", FLANK5, "--flank3", FLANK3,
        "--genome", "/tmp/genome.fsa", "--gene-name", "Fus3",
    ])
    assert "codon_selection_positions" not in base


# --- --orf-plus wiring ---------------------------------------------------

def test_orf_plus_flag_parses_into_kwarg():
    kwargs, _ = cli.build_kwargs([
        "--orf", "o.fa", "--orf-plus", "op.fa", "--gene-name", "G",
        "--genome", "g.fsa", "--output", "out"])
    assert kwargs["orf_plus_file_path"] == "op.fa"


def test_orf_plus_absent_by_default():
    kwargs, _ = cli.build_kwargs([
        "--orf", "o.fa", "--flank5", "a.fa", "--flank3", "b.fa",
        "--gene-name", "G", "--genome", "g.fsa", "--output", "out"])
    assert not cli._is_set(kwargs.get("orf_plus_file_path"))


def test_orf_plus_is_a_per_orf_key():
    assert "orf_plus_file_path" in cli.PER_ORF_KEYS


def test_manifest_accepts_an_orf_plus_column(tmp_path):
    manifest = tmp_path / "orfs.tsv"
    manifest.write_text("gene\torf\torf_plus\nFus3\tfus3.fa\tfus3_plus.fa\n")
    rows = cli._load_manifest(str(manifest))
    assert rows[0]["orf_plus_file_path"].endswith("fus3_plus.fa")


def test_folder_discovery_recognizes_an_orf_plus_suffix(tmp_path):
    (tmp_path / "Fus3_coding.fa").write_text(">x\nATG\n")
    (tmp_path / "Fus3_orfPlus.fa").write_text(">x\nAAATGTTT\n")
    orfs, skipped = cli.discover_orf_folder(str(tmp_path))
    assert skipped == []
    assert len(orfs) == 1
    assert orfs[0]["orf_plus_file_path"].endswith("Fus3_orfPlus.fa")
    assert orfs[0]["orf_file_path"].endswith("Fus3_coding.fa")


def test_orf_plus_suffix_does_not_swallow_the_plain_orf_suffix(tmp_path):
    """'_orf' must not match '_orfPlus', or the ORF and the combined file collide."""
    (tmp_path / "Fus3_orf.fa").write_text(">x\nATG\n")
    (tmp_path / "Fus3_orfPlus.fa").write_text(">x\nAAATGTTT\n")
    orfs, _ = cli.discover_orf_folder(str(tmp_path))
    assert orfs[0]["orf_file_path"].endswith("Fus3_orf.fa")
    assert orfs[0]["orf_plus_file_path"].endswith("Fus3_orfPlus.fa")
