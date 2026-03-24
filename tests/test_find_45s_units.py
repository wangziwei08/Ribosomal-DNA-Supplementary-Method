from pathlib import Path

import find_45s_units as mod


def test_calls_plus_strand_unit_with_58s_overlapping_28s():
    feats = [
        mod.Feature("chr1", 100, 1900, "+", "18S_rRNA"),
        mod.Feature("chr1", 2100, 2253, "+", "5_8S_rRNA"),
        mod.Feature("chr1", 2090, 5900, "+", "28S_rRNA"),
    ]
    units = mod.find_45s_units(feats)
    assert len(units) == 1
    assert units[0].chrom == "chr1"
    assert units[0].strand == "+"
    assert units[0].start == 100
    assert units[0].end == 5900


def test_calls_minus_strand_unit_when_components_share_strand():
    feats = [
        mod.Feature("chr2", 5000, 6800, "-", "18S_rRNA"),
        mod.Feature("chr2", 7100, 7253, "-", "5_8S_rRNA"),
        mod.Feature("chr2", 7090, 10850, "-", "28S_rRNA"),
    ]
    units = mod.find_45s_units(feats)
    assert len(units) == 1
    assert units[0].strand == "-"


def test_rejects_mixed_strands_and_long_distance():
    feats = [
        mod.Feature("chr3", 100, 1900, "+", "18S_rRNA"),
        mod.Feature("chr3", 2100, 2253, "-", "5_8S_rRNA"),
        mod.Feature("chr3", 2090, 5900, "+", "28S_rRNA"),
        mod.Feature("chr3", 20000, 21800, "+", "18S_rRNA"),
        mod.Feature("chr3", 22100, 22253, "+", "5_8S_rRNA"),
        mod.Feature("chr3", 42090, 45900, "+", "28S_rRNA"),
    ]
    units = mod.find_45s_units(feats, max_span=10000)
    assert units == []


def test_summary_text_mentions_28s_as_26s_equivalent(tmp_path: Path):
    out = tmp_path / "summary.txt"
    mod.write_summary(
        out,
        {
            "hap1": {"total_units": 2, "chromosomes": [("chr14_1", 2, 3998, 23934, "-")]},
            "hap2": {"total_units": 3, "chromosomes": [("chr14_2", 3, 7798, 37913, "-")]},
        },
    )
    text = out.read_text(encoding="utf-8")
    assert "28S_rRNA as the 26S-equivalent LSU annotation" in text
    assert "hap1" in text
    assert "chr14_2" in text
