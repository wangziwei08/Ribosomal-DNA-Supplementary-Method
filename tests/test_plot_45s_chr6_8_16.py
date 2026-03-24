from pathlib import Path

import plot_45s_chr6_8_16 as mod


def test_bin_counts_assign_units_by_midpoint():
    units = [(100, 200), (900_000, 900_100), (1_100_000, 1_100_100)]
    bins = mod.build_bins(units, chrom_len=2_000_000, bin_size=1_000_000)
    assert bins == [2, 1]


def test_svg_contains_fixed_y_axis_labels(tmp_path: Path):
    out = tmp_path / "45s.svg"
    mod.render_svg(out, panel_data={})
    text = out.read_text(encoding="utf-8")
    assert ">0<" in text
    assert ">5<" in text
    assert ">10<" in text
    assert "45S unit count per window" not in text


def test_shared_global_scale_shortens_non_max_chromosomes():
    width = 1000
    left = 100
    right = 100
    full = width - left - right
    shorter = mod.scaled_width(15, 20, full)
    assert shorter == full * 0.75
