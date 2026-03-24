from pathlib import Path

import build_5s_98_outputs as mod


def test_build_bins_counts_hits_by_midpoint():
    hits = [(1, 100), (200_001, 200_100), (250_001, 250_100), (510_000, 510_100)]
    bins = mod.build_bins(hits, chrom_len=1_000_000, bin_size=250_000)
    assert bins == [2, 1, 1, 0]


def test_methods_results_text_contains_core_phrasing(tmp_path: Path):
    out = tmp_path / "draft.md"
    mod.write_methods_results(
        out,
        counts={
            ("hap1", "chr2_1"): 591,
            ("hap1", "chr12_1"): 1687,
            ("hap2", "chr2_2"): 1875,
            ("hap2", "chr12_2"): 1275,
        },
    )
    text = out.read_text(encoding="utf-8")
    assert "98% sequence identity" in text
    assert "100% query coverage" in text
    assert "chr2_1" in text
    assert "chr12_2" in text


def test_svg_contains_fixed_density_axis_labels(tmp_path: Path):
    out = tmp_path / "5s.svg"
    mod.render_svg(
        out,
        hits_by_panel={
            "hap1": {"chr2_1": []},
            "hap2": {},
        },
        clusters={"hap1": {"chr2_1": []}, "hap2": {}},
    )
    text = out.read_text(encoding="utf-8")
    assert ">0<" in text
    assert ">100<" in text
    assert ">500<" in text
