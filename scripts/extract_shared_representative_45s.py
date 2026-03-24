from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


BASE = Path("/home/dell/disk2_16T/fojia/Rfam")
TMP = BASE / "tmp"
KEEP = {
    "hap1": {"chr6_1", "chr8_1", "chr16_1"},
    "hap2": {"chr6_2", "chr8_2", "chr16_2"},
}


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def load_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def collect_sequences() -> tuple[dict[str, list[tuple[str, str, int, int]]], dict[str, list[tuple[str, str, int, int]]]]:
    genomes = {
        "hap1": load_fasta(BASE / "fojia.hap1.filtered.fa"),
        "hap2": load_fasta(BASE / "fojia.hap2.filtered.fa"),
    }
    collected = {"hap1": defaultdict(list), "hap2": defaultdict(list)}
    for genome in ("hap1", "hap2"):
        with (TMP / f"{genome}_45S_units.tsv").open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                chrom = row["chromosome"]
                if chrom not in KEEP[genome]:
                    continue
                start = int(row["unit_start"])
                end = int(row["unit_end"])
                strand = row["strand"]
                seq = genomes[genome][chrom][start - 1 : end]
                if strand == "-":
                    seq = revcomp(seq)
                collected[genome][seq.upper()].append((chrom, strand, start, end))
    return collected["hap1"], collected["hap2"]


def pick_shared_main_type(hap1: dict, hap2: dict):
    shared = set(hap1) & set(hap2)
    if not shared:
        raise ValueError("No shared exact 45S sequence found between hap1 and hap2.")
    best = None
    for seq in shared:
        score = (len(hap1[seq]) + len(hap2[seq]), len(hap1[seq]), len(hap2[seq]), len(seq))
        if best is None or score > best[0]:
            best = (score, seq)
    seq = best[1]
    return seq, len(hap1[seq]), len(hap2[seq]), hap1[seq][0], hap2[seq][0]


def write_summary(out: Path, summary: dict) -> None:
    lines = [
        "Shared representative 45S sequence",
        "most abundant exact 45S sequence shared by hap1 and hap2",
        "",
        "length\t" + str(summary["length"]),
        "hap1_count\t" + str(summary["hap1_count"]),
        "hap2_count\t" + str(summary["hap2_count"]),
        "total_count\t" + str(summary["total_count"]),
        "hap1_meta\t" + "\t".join(map(str, summary["hap1_meta"])),
        "hap2_meta\t" + "\t".join(map(str, summary["hap2_meta"])),
    ]
    out.write_text("\n".join(lines), encoding="utf-8")


def run() -> None:
    hap1, hap2 = collect_sequences()
    seq, h1n, h2n, meta1, meta2 = pick_shared_main_type(hap1, hap2)
    header = (
        f"shared_representative_45S|len={len(seq)}|total={h1n+h2n}|hap1={h1n}|hap2={h2n}"
        f"|hap1_locus={meta1[0]}:{meta1[2]}-{meta1[3]}({meta1[1]})"
        f"|hap2_locus={meta2[0]}:{meta2[2]}-{meta2[3]}({meta2[1]})"
    )
    out_fa = TMP / "shared_representative_45S.fa"
    with out_fa.open("w") as fo:
        fo.write(f">{header}\n")
        fo.write(wrap_fasta(seq) + "\n")
    write_summary(
        TMP / "shared_representative_45S_summary.txt",
        {
            "length": len(seq),
            "hap1_count": h1n,
            "hap2_count": h2n,
            "total_count": h1n + h2n,
            "hap1_meta": meta1,
            "hap2_meta": meta2,
        },
    )
    print(out_fa)


if __name__ == "__main__":
    run()
