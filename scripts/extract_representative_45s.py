from __future__ import annotations

import csv
from collections import Counter
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


def write_summary(out: Path, records: list[dict]) -> None:
    lines = ["Representative 45S sequences", "most abundant exact 45S sequence type", ""]
    if records:
        header = "\t".join(records[0].keys())
        lines.append(header)
        for rec in records:
            lines.append("\t".join(str(rec[k]) for k in rec))
    out.write_text("\n".join(lines), encoding="utf-8")


def run() -> None:
    genomes = {
        "hap1": load_fasta(BASE / "fojia.hap1.filtered.fa"),
        "hap2": load_fasta(BASE / "fojia.hap2.filtered.fa"),
    }
    fasta_records: list[tuple[str, str]] = []
    summary_records: list[dict] = []

    for genome in ("hap1", "hap2"):
        counter: Counter[str] = Counter()
        best_meta: dict[str, tuple[str, str, int, int]] = {}
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
                seq = seq.upper()
                counter[seq] += 1
                best_meta.setdefault(seq, (chrom, strand, start, end))

        if not counter:
            continue
        seq, count = counter.most_common(1)[0]
        chrom, strand, start, end = best_meta[seq]
        header = (
            f"{genome}_representative_45S|chrom={chrom}|strand={strand}|start={start}|end={end}"
            f"|len={len(seq)}|copy_count={count}"
        )
        fasta_records.append((header, seq))
        summary_records.append(
            {
                "genome": genome,
                "chromosome": chrom,
                "strand": strand,
                "start": start,
                "end": end,
                "length": len(seq),
                "copy_count": count,
            }
        )

    out_fa = TMP / "representative_45S.fa"
    with out_fa.open("w") as fo:
        for header, seq in fasta_records:
            fo.write(f">{header}\n")
            fo.write(wrap_fasta(seq) + "\n")

    write_summary(TMP / "representative_45S_summary.txt", summary_records)
    print(out_fa)


if __name__ == "__main__":
    run()
