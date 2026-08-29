"""Unit tests for the pure sequence-manipulation helpers (no BLAST required)."""

from pam_scanning import library as L


def test_reverse_complement_basic():
    assert L.reverseComplement("ATGC") == "GCAT"
    assert L.reverseComplement("AAAA") == "TTTT"


def test_reverse_complement_preserves_case():
    assert L.reverseComplement("atgc") == "gcat"


def test_reverse_complement_is_involution():
    seq = "ATGCGTACCGGATCAGT"
    assert L.reverseComplement(L.reverseComplement(seq)) == seq


def test_complement_no_reversal():
    assert L.complement("ATGC") == "TACG"


def test_count_mismatches():
    assert L.countMismatches("AAAA", "AAAA") == 0
    assert L.countMismatches("AAAA", "AAGA") == 1
    assert L.countMismatches("ACGT", "TGCA") == 4


def test_fasta_wraps_at_60():
    seq = "A" * 130
    wrapped = L.fasta(seq)
    lines = [ln for ln in wrapped.split("\n") if ln]
    assert lines[0] == "A" * 60
    assert lines[1] == "A" * 60
    assert lines[2] == "A" * 10
    # No wrapped line exceeds 60 characters.
    assert all(len(ln) <= 60 for ln in lines)


def test_mark_silencers_lowercases_unchanged_bases():
    # Identical sequences -> everything lowercased (no mutations highlighted).
    assert L.markSilencers("ATGC", "ATGC") == "atgc"
    # A single change is kept uppercase to flag the silencing mutation.
    assert L.markSilencers("ATTC", "ATGC") == "atTc"


# --- codon-position parsing (shared by CLI --codon-positions and the GUI picker) ---

from pam_scanning.chimeras import parse_codon_positions as _pcp


def test_parse_codon_positions_singletons_and_ranges():
    assert _pcp("52, 89, 100-105") == [52, 89, 100, 101, 102, 103, 104, 105]


def test_parse_codon_positions_dedupes_and_sorts():
    assert _pcp("3 3 1  2-4") == [1, 2, 3, 4]


def test_parse_codon_positions_clamps_to_length():
    assert _pcp("0, 5, 300", n=10) == [5]


def test_parse_codon_positions_ignores_junk_and_empty():
    assert _pcp("") == []
    assert _pcp(None) == []
    assert _pcp("1-2-3, 7, abc") == [7]   # malformed range dropped, valid kept


# --- ORF + flanks in one file -------------------------------------------

def test_split_orf_plus_returns_the_flanks_around_the_orf():
    from pam_scanning.chimeras import split_orf_plus
    flank5, orf, flank3 = "AAAACCCC", "ATGGGGTAA", "TTTTGGGG"
    assert split_orf_plus(flank5 + orf + flank3, orf) == (flank5, flank3)


def test_split_orf_plus_is_case_insensitive_and_uppercases():
    from pam_scanning.chimeras import split_orf_plus
    assert split_orf_plus("aaaaATGTAAtttt", "atgtaa") == ("AAAA", "TTTT")


def test_split_orf_plus_does_not_assume_a_flank_length():
    """Any flank size works, so 30 bp files are handled as well as 100 bp ones."""
    from pam_scanning.chimeras import split_orf_plus
    orf = "ATG" + "GCT" * 40 + "TAA"
    for n in (5, 30, 100, 250):
        f5, f3 = "A" * n, "T" * n
        assert split_orf_plus(f5 + orf + f3, orf) == (f5, f3)


def test_split_orf_plus_rejects_a_missing_orf():
    import pytest
    from pam_scanning.chimeras import split_orf_plus
    with pytest.raises(ValueError, match="not found"):
        split_orf_plus("AAAACCCCTTTT", "ATGGGGTAA")


def test_split_orf_plus_rejects_an_ambiguous_orf():
    import pytest
    from pam_scanning.chimeras import split_orf_plus
    orf = "ATGGGGTAA"
    with pytest.raises(ValueError, match="occurs 2 times"):
        split_orf_plus("AAAA" + orf + "TTTT" + orf, orf)


def test_split_orf_plus_rejects_a_zero_length_flank():
    import pytest
    from pam_scanning.chimeras import split_orf_plus
    orf = "ATGGGGTAA"
    with pytest.raises(ValueError, match="zero length"):
        split_orf_plus(orf + "TTTT", orf)
    with pytest.raises(ValueError, match="zero length"):
        split_orf_plus("AAAA" + orf, orf)
