from pathlib import Path

import extract_shared_representative_45s as mod


def test_pick_most_abundant_shared_sequence():
    hap1 = {"AAA": [("chr1", "+", 1, 3), ("chr1", "+", 10, 12)], "CCC": [("chr1", "+", 20, 22)]}
    hap2 = {"AAA": [("chr2", "-", 5, 7)], "GGG": [("chr2", "+", 8, 10)]}
    seq, h1n, h2n, meta1, meta2 = mod.pick_shared_main_type(hap1, hap2)
    assert seq == "AAA"
    assert h1n == 2
    assert h2n == 1
    assert meta1 == ("chr1", "+", 1, 3)
    assert meta2 == ("chr2", "-", 5, 7)


def test_summary_mentions_shared_exact_sequence(tmp_path: Path):
    out = tmp_path / "summary.txt"
    mod.write_summary(
        out,
        {
            "length": 5798,
            "hap1_count": 4,
            "hap2_count": 5,
            "total_count": 9,
            "hap1_meta": ("chr16_1", "+", 1, 5798),
            "hap2_meta": ("chr16_2", "+", 2, 5799),
        },
    )
    text = out.read_text(encoding="utf-8")
    assert "most abundant exact 45S sequence shared by hap1 and hap2" in text
    assert "total_count" in text
