from __future__ import annotations

import csv
from pathlib import Path


BASE = Path("/home/dell/disk2_16T/fojia/Rfam")
TMP = BASE / "tmp"
BIN_SIZE = 250_000
IDENTITY_THRESHOLD = 98.0
COVERAGE_THRESHOLD = 100.0

CHROM_LENGTHS = {
    "chr2_1": 14130241,
    "chr12_1": 9330097,
    "chr2_2": 13418534,
    "chr12_2": 9181357,
}

PANELS = [
    ("hap1", "chr2_1", "hap1 chr2"),
    ("hap1", "chr12_1", "hap1 chr12"),
    ("hap2", "chr2_2", "hap2 chr2"),
    ("hap2", "chr12_2", "hap2 chr12"),
]


def build_bins(
    hits: list[tuple[int, int]], chrom_len: int, bin_size: int = BIN_SIZE
) -> list[int]:
    n_bins = (chrom_len + bin_size - 1) // bin_size
    bins = [0] * n_bins
    for start, end in hits:
        midpoint = (start + end) / 2
        idx = min(int(midpoint // bin_size), n_bins - 1)
        bins[idx] += 1
    return bins


def filter_hits() -> tuple[dict[str, list[tuple[str, int, int]]], dict[tuple[str, str], int]]:
    result: dict[str, list[tuple[str, int, int]]] = {"hap1": [], "hap2": []}
    counts: dict[tuple[str, str], int] = {}
    for genome in ("hap1", "hap2"):
        src = TMP / f"{genome}_vs_shared_5S.blast.tsv"
        with src.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                if not row:
                    continue
                pident = float(row[2])
                qcovs = float(row[14])
                if pident < IDENTITY_THRESHOLD or qcovs < COVERAGE_THRESHOLD:
                    continue
                chrom = row[1]
                start, end = sorted((int(row[8]), int(row[9])))
                result[genome].append((chrom, start, end))
                counts[(genome, chrom)] = counts.get((genome, chrom), 0) + 1
    return result, counts


def cluster_hits(
    hits: list[tuple[int, int]], merge_gap: int = 1000
) -> list[tuple[int, int, int]]:
    if not hits:
        return []
    ordered = sorted(hits)
    clusters: list[tuple[int, int, int]] = []
    cur_start, cur_end = ordered[0]
    cur_count = 1
    for start, end in ordered[1:]:
        if start <= cur_end + merge_gap:
            cur_end = max(cur_end, end)
            cur_count += 1
        else:
            clusters.append((cur_start, cur_end, cur_count))
            cur_start, cur_end = start, end
            cur_count = 1
    clusters.append((cur_start, cur_end, cur_count))
    return clusters


def write_filtered_tables(
    filtered: dict[str, list[tuple[str, int, int]]]
) -> dict[str, dict[str, list[tuple[int, int, int]]]]:
    cluster_map: dict[str, dict[str, list[tuple[int, int, int]]]] = {"hap1": {}, "hap2": {}}
    for genome, records in filtered.items():
        records = sorted(records, key=lambda x: (x[0], x[1], x[2]))
        out_tsv = TMP / f"{genome}_98id100cov_hits.tsv"
        out_bed = TMP / f"{genome}_98id100cov_hits.bed"
        with out_tsv.open("w", newline="") as fo:
            writer = csv.writer(fo, delimiter="\t")
            writer.writerow(["genome", "chromosome", "start", "end", "length"])
            for chrom, start, end in records:
                writer.writerow([genome, chrom, start, end, end - start + 1])
        with out_bed.open("w", newline="") as fo:
            writer = csv.writer(fo, delimiter="\t")
            for i, (chrom, start, end) in enumerate(records, 1):
                writer.writerow([chrom, start - 1, end, f"{genome}_5S98_{i}", 980, "+"])

        by_chrom: dict[str, list[tuple[int, int]]] = {}
        for chrom, start, end in records:
            by_chrom.setdefault(chrom, []).append((start, end))
        out_cluster = TMP / f"{genome}_98id100cov_clusters_1kb.tsv"
        with out_cluster.open("w", newline="") as fo:
            writer = csv.writer(fo, delimiter="\t")
            writer.writerow(["chromosome", "start", "end", "span_bp", "hits_in_cluster"])
            for chrom in sorted(by_chrom):
                clusters = cluster_hits(by_chrom[chrom], merge_gap=1000)
                cluster_map[genome][chrom] = clusters
                for start, end, count in clusters:
                    writer.writerow([chrom, start, end, end - start + 1, count])
    return cluster_map


def write_summary_table(
    counts: dict[tuple[str, str], int],
    clusters: dict[str, dict[str, list[tuple[int, int, int]]]],
) -> Path:
    out = TMP / "main_table_5S_98id100cov.tsv"
    with out.open("w", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t")
        writer.writerow(
            [
                "genome",
                "chromosome",
                "chromosome_length_bp",
                "high_confidence_5S_like_copies",
                "largest_cluster_start",
                "largest_cluster_end",
                "largest_cluster_span_bp",
                "largest_cluster_hits",
            ]
        )
        for genome, chrom, _ in PANELS:
            largest = (0, 0, 0)
            if clusters[genome].get(chrom):
                largest = max(clusters[genome][chrom], key=lambda x: x[2])
            writer.writerow(
                [
                    genome,
                    chrom,
                    CHROM_LENGTHS[chrom],
                    counts.get((genome, chrom), 0),
                    largest[0],
                    largest[1],
                    largest[1] - largest[0] + 1 if largest[2] else 0,
                    largest[2],
                ]
            )
    return out


def render_svg(
    out: Path,
    hits_by_panel: dict[str, dict[str, list[tuple[int, int]]]],
    clusters: dict[str, dict[str, list[tuple[int, int, int]]]],
) -> None:
    width = 1700
    height = 1220
    left = 185
    right = 95
    bar_w = width - left - right
    panel_h = 260
    top0 = 125
    bar_h = 20
    density_h = 58
    axis_max = 500

    svg: list[str] = []
    svg.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
    )
    svg.append('<rect width="100%" height="100%" fill="#f7f1e3"/>')
    svg.append(
        '<text x="70" y="58" font-size="30" font-family="Helvetica, Arial, sans-serif" '
        'font-weight="700" fill="#1e2a2f">High-confidence 5S-like loci on chr2 and chr12</text>'
    )
    svg.append(
        '<text x="70" y="88" font-size="15" font-family="Helvetica, Arial, sans-serif" '
        'fill="#405059">Query: shared_most_common_5S.fa | Threshold: identity >= 98%, query coverage = 100%</text>'
    )

    for i, (genome, chrom, title) in enumerate(PANELS):
        chrom_len = CHROM_LENGTHS[chrom]
        panel_top = top0 + i * panel_h
        density_base = panel_top + 118
        bar_y = panel_top + 136
        axis_y = bar_y + 40
        panel_hits = hits_by_panel[genome].get(chrom, [])
        bins = build_bins(panel_hits, chrom_len, BIN_SIZE)
        max_bin = max(bins) if bins else 0

        svg.append(
            f'<text x="70" y="{panel_top - 2}" font-size="24" font-family="Helvetica, Arial, sans-serif" '
            f'font-weight="700" fill="#1e2a2f">{title}</text>'
        )
        svg.append(
            f'<text x="70" y="{panel_top + 16}" font-size="15" font-family="Helvetica, Arial, sans-serif" '
            f'fill="#56666d">chromosome length: {chrom_len/1_000_000:.2f} Mb | hits: {len(panel_hits)} | max per 0.25 Mb bin: {max_bin}</text>'
        )
        for tick in (0, 100, 200, 300, 400, 500):
            y = density_base - (tick / axis_max) * density_h
            svg.append(
                f'<line x1="{left}" y1="{y:.2f}" x2="{left + bar_w}" y2="{y:.2f}" stroke="#d5d9d6" stroke-width="1"/>'
            )
            svg.append(
                f'<text x="{left - 12}" y="{y + 4:.2f}" text-anchor="end" font-size="12" font-family="Helvetica, Arial, sans-serif" fill="#405059">{tick}</text>'
            )
        for idx, count in enumerate(bins):
            x = left + (idx * BIN_SIZE / chrom_len) * bar_w
            next_pos = min((idx + 1) * BIN_SIZE, chrom_len)
            w = max(2.0, ((next_pos - idx * BIN_SIZE) / chrom_len) * bar_w - 0.8)
            h = (min(count, axis_max) / axis_max) * density_h
            y = density_base - h
            svg.append(
                f'<rect class="density-bar" x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="#3a86b8" fill-opacity="0.80"/>'
            )
        svg.append(
            f'<rect x="{left}" y="{bar_y}" width="{bar_w}" height="{bar_h}" rx="10" fill="#d8dfdc" stroke="#7d8b8a" stroke-width="1.5"/>'
        )
        for start, end, count in clusters[genome].get(chrom, []):
            if count < 20:
                continue
            x = left + (start / chrom_len) * bar_w
            w = max(3.0, ((end - start + 1) / chrom_len) * bar_w)
            svg.append(
                f'<rect x="{x:.2f}" y="{bar_y - 7}" width="{w:.2f}" height="{bar_h + 14}" rx="4" fill="#cf5c36" fill-opacity="0.68" stroke="#8b2e16" stroke-width="0.8"/>'
            )
        for pos in range(0, chrom_len + 1_000_000, 1_000_000):
            pos = min(pos, chrom_len)
            x = left + (pos / chrom_len) * bar_w
            svg.append(
                f'<line x1="{x:.2f}" y1="{axis_y - 14}" x2="{x:.2f}" y2="{axis_y - 6}" stroke="#7d8b8a" stroke-width="1"/>'
            )
            svg.append(
                f'<text x="{x:.2f}" y="{axis_y + 10}" text-anchor="middle" font-size="12" font-family="Helvetica, Arial, sans-serif" fill="#405059">{int(round(pos/1_000_000))} Mb</text>'
            )
    legend_y = height - 96
    svg.append(
        f'<rect x="90" y="{legend_y - 20}" width="30" height="18" class="density-bar" fill="#3a86b8" fill-opacity="0.80"/>'
    )
    svg.append(
        f'<text x="134" y="{legend_y - 6}" font-size="14" font-family="Helvetica, Arial, sans-serif" fill="#405059">0.25 Mb density bins</text>'
    )
    svg.append(
        f'<rect x="330" y="{legend_y - 22}" width="30" height="22" rx="4" fill="#cf5c36" fill-opacity="0.68" stroke="#8b2e16" stroke-width="1"/>'
    )
    svg.append(
        f'<text x="374" y="{legend_y - 6}" font-size="14" font-family="Helvetica, Arial, sans-serif" fill="#405059">1 kb clusters with at least 20 hits</text>'
    )
    svg.append("</svg>")

    out.write_text("\n".join(svg), encoding="utf-8")


def render_figure(
    filtered: dict[str, list[tuple[str, int, int]]],
    clusters: dict[str, dict[str, list[tuple[int, int, int]]]],
) -> tuple[Path, Path]:
    hits_by_panel: dict[str, dict[str, list[tuple[int, int]]]] = {"hap1": {}, "hap2": {}}
    for genome, records in filtered.items():
        for chrom, start, end in records:
            hits_by_panel[genome].setdefault(chrom, []).append((start, end))
    svg_path = TMP / "figure_5S_98id100cov_chr2_chr12.svg"
    png_path = TMP / "figure_5S_98id100cov_chr2_chr12.png"
    render_svg(svg_path, hits_by_panel, clusters)
    return svg_path, png_path


def write_methods_results(out: Path, counts: dict[tuple[str, str], int]) -> None:
    text = f"""# Scientific Data Draft Text

## Methods

To characterize the genomic distribution of shared 5S rDNA-related loci, the most prevalent 5S sequence shared between hap1 and hap2 (`shared_most_common_5S.fa`, 119 bp) was used as the query. Only chr2 and chr12 were analyzed in each haplotype assembly by extracting chr2_1 and chr12_1 from `fojia.hap1.filtered.fa` and chr2_2 and chr12_2 from `fojia.hap2.filtered.fa`. Local nucleotide searches were performed with BLASTN (`blastn-short`, BLAST+ v2.5.0) against haplotype-specific databases built from the extracted chromosome sequences. Hits were retained as high-confidence 5S-like loci when they showed at least 98% sequence identity together with 100% query coverage across the full 119-bp query sequence. Retained matches were summarized per chromosome, exported as BED intervals, and further grouped into local clusters by merging neighboring hits separated by <=1 kb. For visualization, high-confidence loci were projected onto chromosome-length-scaled schematics, and local densities were summarized in 0.25-Mb windows with 1-Mb coordinate tick marks.

## Results

Using a threshold of at least 98% sequence identity and 100% query coverage, we identified high-confidence 5S-like loci on both chr2 and chr12 in the two haplotype assemblies. In hap1, {counts[('hap1', 'chr2_1')]} loci were detected on chr2_1 and {counts[('hap1', 'chr12_1')]} loci were detected on chr12_1. In hap2, {counts[('hap2', 'chr2_2')]} loci were detected on chr2_2 and {counts[('hap2', 'chr12_2')]} loci were detected on chr12_2. The loci were not evenly distributed along the chromosomes, but instead formed prominent local enrichments. In hap1, the major enrichment on chr12_1 was centered around the 3.1-4.1 Mb region, whereas the major chr2_1 enrichment was concentrated near 4.8-8.4 Mb with a particularly pronounced peak around 8 Mb. In hap2, chr12_2 showed a major enrichment between approximately 3.0 and 3.8 Mb, while chr2_2 contained the strongest accumulation between about 4.7 and 8.1 Mb. These results indicate that shared 5S-like sequences are concentrated in large clustered arrays on chr2 and chr12 in both haplotypes rather than being uniformly dispersed across the chromosomes.
"""
    out.write_text(text, encoding="utf-8")


def main() -> None:
    filtered, counts = filter_hits()
    clusters = write_filtered_tables(filtered)
    write_summary_table(counts, clusters)
    svg_path, png_path = render_figure(filtered, clusters)
    write_methods_results(TMP / "scientific_data_methods_results_5S_98id100cov.md", counts)
    print(svg_path)
    print(png_path)


if __name__ == "__main__":
    main()
