from __future__ import annotations

import csv
from pathlib import Path


BASE = Path("/home/dell/disk2_16T/fojia/Rfam")
TMP = BASE / "tmp"
BIN_SIZE = 100_000

PANELS = [
    ("hap1", "chr6_1", "Hap1 chr6", 17739215),
    ("hap1", "chr8_1", "Hap1 chr8", 17269009),
    ("hap1", "chr16_1", "Hap1 chr16", 17664751),
    ("hap2", "chr6_2", "Hap2 chr6", 17642119),
    ("hap2", "chr8_2", "Hap2 chr8", 17161392),
    ("hap2", "chr16_2", "Hap2 chr16", 17569814),
]


def scaled_width(chrom_len_mb: float, max_len_mb: float, full_width: float) -> float:
    return full_width * chrom_len_mb / max_len_mb


def load_units() -> dict[tuple[str, str], list[tuple[int, int]]]:
    data: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for genome in ("hap1", "hap2"):
        path = TMP / f"{genome}_45S_units.tsv"
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                key = (genome, row["chromosome"])
                data.setdefault(key, []).append((int(row["unit_start"]), int(row["unit_end"])))
    return data


def build_bins(units: list[tuple[int, int]], chrom_len: int, bin_size: int = BIN_SIZE) -> list[int]:
    n_bins = (chrom_len + bin_size - 1) // bin_size
    bins = [0] * n_bins
    for start, end in units:
        midpoint = (start + end) / 2
        idx = min(int(midpoint // bin_size), n_bins - 1)
        bins[idx] += 1
    return bins


def render_svg(out: Path, panel_data: dict[tuple[str, str], list[tuple[int, int]]] | None = None) -> None:
    panel_data = panel_data or load_units()
    width = 1750
    height = 1520
    left = 190
    right = 95
    track_w = width - left - right
    panel_h = 220
    top0 = 120
    axis_max = 10
    density_h = 68
    chrom_bar_h = 18
    max_len = max(chrom_len for _, _, _, chrom_len in PANELS)

    svg: list[str] = []
    svg.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
    )
    svg.append('<rect width="100%" height="100%" fill="#ffffff"/>')
    svg.append(
        '<text x="70" y="55" font-size="28" font-family="Helvetica, Arial, sans-serif" font-weight="700" fill="#1e2a2f">45S (18S-5.8S-26S) distribution on chr6, chr8 and chr16</text>'
    )

    for i, (genome, chrom, title, chrom_len) in enumerate(PANELS):
        panel_top = top0 + i * panel_h
        density_base = panel_top + 120
        chrom_y = panel_top + 136
        axis_y = chrom_y + 36
        actual_w = scaled_width(chrom_len, max_len, track_w)
        units = panel_data.get((genome, chrom), [])
        bins = build_bins(units, chrom_len)

        svg.append(
            f'<text x="70" y="{panel_top + 4}" font-size="23" font-family="Helvetica, Arial, sans-serif" font-weight="700" fill="#1e2a2f">{title}</text>'
        )
        svg.append(
            f'<text x="70" y="{panel_top + 24}" font-size="16" font-family="Helvetica, Arial, sans-serif" fill="#56666d">chromosome length: {chrom_len/1_000_000:.2f} Mb | complete 45S units: {len(units)}</text>'
        )

        for tick in (0, 5, 10):
            y = density_base - (tick / axis_max) * density_h
            svg.append(
                f'<line x1="{left}" y1="{y:.2f}" x2="{left + track_w}" y2="{y:.2f}" stroke="#d5d9d6" stroke-width="1"/>'
            )
            svg.append(
                f'<text x="{left - 12}" y="{y + 4:.2f}" text-anchor="end" font-size="12" font-family="Helvetica, Arial, sans-serif" fill="#405059">{tick}</text>'
            )

        for idx, count in enumerate(bins):
            x = left + (idx * BIN_SIZE / chrom_len) * actual_w
            next_pos = min((idx + 1) * BIN_SIZE, chrom_len)
            w = max(1.6, ((next_pos - idx * BIN_SIZE) / chrom_len) * actual_w - 0.5)
            clipped = min(count, axis_max)
            h = (clipped / axis_max) * density_h
            svg.append(
                f'<rect x="{x:.2f}" y="{density_base - h:.2f}" width="{w:.2f}" height="{h:.2f}" fill="#3a86b8" fill-opacity="0.82"/>'
            )

        svg.append(
            f'<rect x="{left}" y="{chrom_y}" width="{actual_w:.2f}" height="{chrom_bar_h}" rx="9" fill="#d8dfdc" stroke="#7d8b8a" stroke-width="1.2"/>'
        )

        for start, end in units:
            x = left + (((start + end) / 2) / chrom_len) * actual_w
            svg.append(
                f'<line x1="{x:.2f}" y1="{chrom_y - 7}" x2="{x:.2f}" y2="{chrom_y + chrom_bar_h + 7}" stroke="#cf5c36" stroke-width="1.2" stroke-opacity="0.9"/>'
            )

        for pos in range(0, chrom_len + 1_000_000, 1_000_000):
            pos = min(pos, chrom_len)
            x = left + (pos / max_len) * track_w
            svg.append(
                f'<line x1="{x:.2f}" y1="{axis_y - 12}" x2="{x:.2f}" y2="{axis_y - 5}" stroke="#7d8b8a" stroke-width="1"/>'
            )
            svg.append(
                f'<text x="{x:.2f}" y="{axis_y + 10}" text-anchor="middle" font-size="11.5" font-family="Helvetica, Arial, sans-serif" fill="#405059">{int(round(pos/1_000_000))} Mb</text>'
            )

        svg.append(
            f'<text x="{left + actual_w + 8:.2f}" y="{chrom_y + 13}" font-size="11.5" font-family="Helvetica, Arial, sans-serif" fill="#405059">{chrom_len/1_000_000:.2f} Mb</text>'
        )

    legend_y = height - 52
    svg.append(
        f'<rect x="85" y="{legend_y - 18}" width="28" height="18" fill="#3a86b8" fill-opacity="0.82"/>'
    )
    svg.append(
        f'<text x="125" y="{legend_y - 4}" font-size="19" font-family="Helvetica, Arial, sans-serif" fill="#405059">0.1 Mb density bins</text>'
    )
    svg.append(
        f'<line x1="300" y1="{legend_y - 18}" x2="300" y2="{legend_y}" stroke="#cf5c36" stroke-width="1.2"/>'
    )
    svg.append(
        f'<text x="315" y="{legend_y - 4}" font-size="19" font-family="Helvetica, Arial, sans-serif" fill="#405059">complete 45S units</text>'
    )
    svg.append("</svg>")
    out.write_text("\n".join(svg), encoding="utf-8")


def main() -> None:
    svg_path = TMP / "figure_45S_chr16_chr8_chr6.svg"
    render_svg(svg_path)
    print(svg_path)


if __name__ == "__main__":
    main()
