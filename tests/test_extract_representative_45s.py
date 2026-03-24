from pathlib import Path

import extract_representative_45s as mod


def test_reverse_complement_handles_minus_strand():
    assert mod.revcomp("ACGTN") == "NACGT"


def test_wrap_fasta_breaks_lines():
    wrapped = mod.wrap_fasta("A" * 130, width=60)
    assert wrapped.splitlines() == ["A" * 60, "A" * 60, "A" * 10]


def test_summary_mentions_exact_sequence_type(tmp_path: Path):
    out = tmp_path / "summary.txt"
    mod.write_summary(
        out,
        [
            {
                "genome": "hap1",
                "chromosome": "chr16_1",
                "strand": "+",
                "start": 1,
                "end": 5798,
                "length": 5798,
                "copy_count": 63,
            }
        ],
    )
    text = out.read_text(encoding="utf-8")
    assert "most abundant exact 45S sequence type" in text
    assert "copy_count" in text
