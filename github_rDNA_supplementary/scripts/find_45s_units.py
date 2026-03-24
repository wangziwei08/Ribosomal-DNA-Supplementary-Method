from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


BASE = Path("/home/dell/disk2_16T/fojia/Rfam")
TMP = BASE / "tmp"
TARGET_NAMES = {"18S_rRNA", "5_8S_rRNA", "28S_rRNA"}


@dataclass(frozen=True)
class Feature:
    chrom: str
    start: int
    end: int
    strand: str
    name: str


@dataclass(frozen=True)
class Unit45S:
    chrom: str
    strand: str
    start: int
    end: int
    r18_start: int
    r18_end: int
    r58_start: int
    r58_end: int
    r28_start: int
    r28_end: int


def parse_gff(path: Path) -> list[Feature]:
    feats: list[Feature] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            if "Name=" not in attrs:
                continue
            name = attrs.split("Name=")[1].split(";")[0]
            if name not in TARGET_NAMES:
                continue
            feats.append(
                Feature(
                    chrom=parts[0],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    strand=parts[6],
                    name=name,
                )
            )
    feats.sort(key=lambda x: (x.chrom, x.strand, x.start, x.end, x.name))
    return feats


def _compatible_triplet(
    r18: Feature, r58: Feature, r28: Feature, max_span: int
) -> bool:
    if len({r18.chrom, r58.chrom, r28.chrom}) != 1:
        return False
    if len({r18.strand, r58.strand, r28.strand}) != 1:
        return False
    left = min(r18.start, r58.start, r28.start)
    right = max(r18.end, r58.end, r28.end)
    if right - left + 1 > max_span:
        return False
    if not (r18.end < r28.start):
        return False
    if not (r18.end < r58.start <= r28.end):
        return False
    if not (r28.start <= r58.end <= r28.end):
        return False
    return True


def find_45s_units(features: Iterable[Feature], max_span: int = 10000) -> list[Unit45S]:
    grouped: dict[tuple[str, str], list[Feature]] = {}
    for feat in features:
        grouped.setdefault((feat.chrom, feat.strand), []).append(feat)

    units: list[Unit45S] = []
    for (chrom, strand), feats in grouped.items():
        r18s = [x for x in feats if x.name == "18S_rRNA"]
        r58s = [x for x in feats if x.name == "5_8S_rRNA"]
        r28s = [x for x in feats if x.name == "28S_rRNA"]
        used18: set[int] = set()
        used58: set[int] = set()
        used28: set[int] = set()

        for i18, r18 in enumerate(r18s):
            best: tuple[int, int, int, int] | None = None
            for i58, r58 in enumerate(r58s):
                if i58 in used58:
                    continue
                for i28, r28 in enumerate(r28s):
                    if i28 in used28:
                        continue
                    if not _compatible_triplet(r18, r58, r28, max_span=max_span):
                        continue
                    span = max(r18.end, r58.end, r28.end) - min(r18.start, r58.start, r28.start) + 1
                    score = (span, abs(r58.start - r28.start), abs(r18.end - r58.start), i58, i28)
                    if best is None or score < best:
                        best = score
            if best is None:
                continue
            _, _, _, i58, i28 = best
            if i18 in used18 or i58 in used58 or i28 in used28:
                continue
            r58 = r58s[i58]
            r28 = r28s[i28]
            used18.add(i18)
            used58.add(i58)
            used28.add(i28)
            units.append(
                Unit45S(
                    chrom=chrom,
                    strand=strand,
                    start=min(r18.start, r58.start, r28.start),
                    end=max(r18.end, r58.end, r28.end),
                    r18_start=r18.start,
                    r18_end=r18.end,
                    r58_start=r58.start,
                    r58_end=r58.end,
                    r28_start=r28.start,
                    r28_end=r28.end,
                )
            )
    units.sort(key=lambda x: (x.chrom, x.start, x.end, x.strand))
    return units


def summarize_units(units: list[Unit45S]) -> list[tuple[str, int, int, int, str]]:
    out: list[tuple[str, int, int, int, str]] = []
    by_chrom: dict[str, list[Unit45S]] = {}
    for unit in units:
        by_chrom.setdefault(unit.chrom, []).append(unit)
    for chrom in sorted(by_chrom):
        arr = sorted(by_chrom[chrom], key=lambda x: x.start)
        strands = "".join(sorted(set(x.strand for x in arr)))
        out.append((chrom, len(arr), arr[0].start, arr[-1].end, strands))
    return out


def write_summary(out: Path, summaries: dict[str, dict]) -> None:
    lines = [
        "45S (18S-5.8S-26S) rDNA summary",
        "Operational definition: same-chromosome, same-strand 18S_rRNA + 5_8S_rRNA + 28S_rRNA triplets.",
        "Note: barrnap reports 28S_rRNA as the 26S-equivalent LSU annotation in these eukaryotic assemblies.",
        "The 5.8S call is allowed to overlap the annotated 28S boundary, matching the barrnap annotation pattern.",
        "",
    ]
    for genome in ("hap1", "hap2"):
        info = summaries[genome]
        lines.append(f"{genome}\ttotal_45S_units\t{info['total_units']}")
        for chrom, count, start, end, strands in info["chromosomes"]:
            lines.append(f"{genome}\t{chrom}\t{count}\t{start}\t{end}\t{strands}")
        lines.append("")
    out.write_text("\n".join(lines), encoding="utf-8")


def run() -> None:
    summaries: dict[str, dict] = {}
    for genome, gff in (("hap1", BASE / "hap1.rRNA.gff"), ("hap2", BASE / "hap2.rRNA.gff")):
        feats = parse_gff(gff)
        units = find_45s_units(feats)
        detail_path = TMP / f"{genome}_45S_units.tsv"
        with detail_path.open("w", newline="") as fo:
            w = csv.writer(fo, delimiter="\t")
            w.writerow(
                [
                    "chromosome",
                    "strand",
                    "unit_start",
                    "unit_end",
                    "unit_span_bp",
                    "18S_start",
                    "18S_end",
                    "5.8S_start",
                    "5.8S_end",
                    "26S_equiv_start",
                    "26S_equiv_end",
                ]
            )
            for u in units:
                w.writerow(
                    [
                        u.chrom,
                        u.strand,
                        u.start,
                        u.end,
                        u.end - u.start + 1,
                        u.r18_start,
                        u.r18_end,
                        u.r58_start,
                        u.r58_end,
                        u.r28_start,
                        u.r28_end,
                    ]
                )
        chrom_summary = summarize_units(units)
        sum_path = TMP / f"{genome}_45S_summary.tsv"
        with sum_path.open("w", newline="") as fo:
            w = csv.writer(fo, delimiter="\t")
            w.writerow(["chromosome", "unit_count", "first_unit_start", "last_unit_end", "strands"])
            w.writerows(chrom_summary)
        summaries[genome] = {
            "total_units": len(units),
            "chromosomes": chrom_summary,
        }

    write_summary(TMP / "45S_summary.txt", summaries)
    print(TMP / "45S_summary.txt")


if __name__ == "__main__":
    run()
