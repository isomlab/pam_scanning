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


# --- homology arms and short flanks -------------------------------------
# A 5' flank shorter than the arm used to drive the left slice index negative,
# so Python read from the END of the sequence and the arm came back EMPTY. The
# primer was then nothing but its suffix. These pin the clamped behaviour.

def _assembled(flank5, orf_len=300, flank3=100):
    orf = "ATG" + "GCT" * ((orf_len - 6) // 3) + "TAA"
    return "A" * flank5 + orf + "T" * flank3, flank5


def test_full_arms_when_the_flank_is_long_enough():
    from pam_scanning.chimeras import homology_arms
    seq, i = _assembled(100)
    left, right = homology_arms(seq, i, 70, 70)
    assert len(left) == 70 and len(right) == 70


def test_short_flank_shortens_the_arm_instead_of_emptying_it():
    from pam_scanning.chimeras import homology_arms
    seq, i = _assembled(30)
    left, _right = homology_arms(seq, i, 70, 70)
    assert left, "a short 5' flank must not produce an empty homology arm"
    assert len(left) == 33          # the 30 nt flank plus the 3 nt of the codon


def test_arm_never_wraps_to_the_end_of_the_sequence():
    """The old negative index read the tail of the sequence, which is worse
    than a short arm: it is the wrong sequence entirely."""
    from pam_scanning.chimeras import homology_arms
    seq, i = _assembled(30)
    left, _ = homology_arms(seq, i, 70, 70)
    assert set(left) <= set("AATGC"), left
    assert not left.endswith("T" * 10)


def test_arms_shorten_smoothly_as_the_flank_shrinks():
    from pam_scanning.chimeras import homology_arms
    lengths = []
    for flank5 in (100, 70, 50, 30, 10, 0):
        seq, i = _assembled(flank5)
        left, _ = homology_arms(seq, i, 70, 70)
        lengths.append(len(left))
    assert lengths == sorted(lengths, reverse=True)
    assert lengths[0] == 70 and lengths[-1] == 3


def test_right_arm_clamps_at_the_end_too():
    from pam_scanning.chimeras import homology_arms
    seq, i = _assembled(100, flank3=10)
    insertion = len(seq) - 10 - 3   # near the 3' end
    _left, right = homology_arms(seq, insertion, 70, 70)
    assert 0 < len(right) <= 70
