#!/usr/bin/env python
import argparse
import pysam
import re
from pathlib import Path
from itertools import chain
from dataclasses import dataclass
from typing import Any


@dataclass(slots=True)
class PairedData:
    rnalib: Any
    genome: Any


def fancy_samtools_string_key(string, _re=re.compile("([0-9]+)")):
    return tuple(int(part) if part.isdigit() else part for part in _re.split(string))


def batch_reads(inbam: Path):
    inbam = pysam.AlignmentFile(inbam, "rb")
    iterator = iter(inbam)

    cache = [next(iterator)]
    for read in iterator:
        if read.query_name == cache[-1].query_name:
            cache.append(read)
            continue
        yield cache
        cache = [read]
    yield cache

    inbam.close()


def nms(
    rnalib: list[pysam.AlignedSegment],
    genome: list[pysam.AlignedSegment],
    saveto: pysam.AlignmentFile,
    summary: dict[str, int],
):
    # Keep only top alignments
    asthr = max([x.get_tag("AS") for x in chain(rnalib, genome)], default=0)
    rnalib = [x for x in rnalib if x.get_tag("AS") == asthr]
    genome = [x for x in genome if x.get_tag("AS") == asthr]

    # Run NMS
    match (len(rnalib), len(genome)):
        case (0, 0):  # Unmapped
            raise NotImplementedError("Unreacheable")
        case (1, 1):  # Unique: genome
            summary["Unique: genome"] += 1
            to_write = genome
        # rnalib -----------------------------------------------------
        case (1, 0):  # Unique: rnalib
            summary["Unique: RNA library"] += 1
            to_write = rnalib
        case (n, 0):  # Multimap: rnalib
            assert n > 1
            summary["Multimap: RNA library"] += 1
            to_write = rnalib
        case (1, n):  # NMS: rnalib
            assert n > 1
            summary["NMS: RNA library"] += 1
            to_write = rnalib
        # genome -----------------------------------------------------
        case (0, 1):  # Unique: genome
            summary["Unique: genome"] += 1
            to_write = genome
        case (0, n):  # Multimap: genome
            assert n > 1
            summary["Multimap: genome"] += 1
            to_write = genome
        case (n, 1):  # NMS: genome
            assert n > 1
            summary["NMS: genome"] += 1
            to_write = genome
        # ------------------------------------------------------------

        case (_, _):  # Ambiguous
            summary["NMS: ambiguous"] += 1
            to_write = rnalib

        case _:
            raise NotImplementedError(
                f"Unknown alignments combination: RNA library = {len(rnalib)}, genome = {len(genome)}"
            )

    for record in to_write:
        # We must create a new record, because previous pysam.AlignedSegment is tightly matched with old header
        saveto.write(
            pysam.AlignedSegment.from_dict(record.to_dict(), saveto.header)
        )


def derive_header(rnalib: Path, genome: Path) -> dict:
    headers = []
    for bam in rnalib, genome:
        with pysam.AlignmentFile(bam, "rb") as stream:
            headers.append(stream.header)

    header = {
        "RG": headers[0]["RG"],
        "HD": headers[0]["HD"],
        "PG": headers[0]["PG"],
        "CO": [],
        "SQ": [],
    }
    for h in headers:
        assert h["RG"] == header["RG"] and h["HD"] == header["HD"]

        header["CO"].extend(h["CO"])
        # header["PG"].extend(h["PG"])
        header["SQ"].extend(h["SQ"])
    return header


parser = argparse.ArgumentParser(
    description="Suppress suboptimal iCLIP alignments to genomic/rnalib libraries"
)
parser.add_argument("threads", type=int, help="Number of threads to use")
parser.add_argument("reads", type=Path, help="Sequenced reads (.fastq.gz)")
parser.add_argument("rnalib", type=Path, help="Alignment to the rnalib library (.bam)")
parser.add_argument("genome", type=Path, help="Alignment to the genome (.bam)")
parser.add_argument(
    "logs", type=Path, help="Save NMS summary stats to the given file (.yaml)"
)
parser.add_argument(
    "alignments", type=Path, help="Save NMS results to the given file (.bam)"
)

# args = ["64", "reads.fastq.gz", "rnalib.bam", "genome.bam", "stats.tsv", "nms.bam"]
# parser = parser.parse_args(args)

parser = parser.parse_args()

# Sort by name input files
print("Sorting input reads by name....")

sorted_by_name = PairedData(
    rnalib=parser.rnalib.with_suffix(".sort-names.bam"),
    genome=parser.genome.with_suffix(".sort-names.bam"),
)
for bam, saveto in (parser.rnalib, sorted_by_name.rnalib), (
    parser.genome,
    sorted_by_name.genome,
):
    pysam.sort("-@", str(parser.threads), "-o", saveto.as_posix(), "-n", bam.as_posix())

# Batch reads & run NMS & save results
saveto = pysam.AlignmentFile(
    parser.alignments,
    "wb",
    header=derive_header(sorted_by_name.rnalib, sorted_by_name.genome),
)
summary = {
    "Unique: RNA library": 0,
    "Multimap: RNA library": 0,
    "NMS: RNA library": 0,
    "Unique: genome": 0,
    "Multimap: genome": 0,
    "NMS: genome": 0,
    "NMS: ambiguous": 0,
}

batches = PairedData(
    rnalib=batch_reads(sorted_by_name.rnalib), genome=batch_reads(sorted_by_name.genome)
)

cursor = PairedData(rnalib=next(batches.rnalib, None), genome=next(batches.genome, None))
while True:
    if cursor.rnalib is None or cursor.genome is None:
        break

    lib_key = fancy_samtools_string_key(cursor.rnalib[0].query_name)
    gn_key = fancy_samtools_string_key(cursor.genome[0].query_name)

    if lib_key == gn_key:
        nms(cursor.rnalib, cursor.genome, saveto, summary)
        cursor.rnalib = next(batches.rnalib, None)
        cursor.genome = next(batches.genome, None)
    elif lib_key < gn_key:
        nms(cursor.rnalib, [], saveto, summary)
        cursor.rnalib = next(batches.rnalib, None)
    else:
        nms([], cursor.genome, saveto, summary)
        cursor.genome = next(batches.genome, None)

# Left-over alignments
while cursor.rnalib is not None:
    nms(cursor.rnalib, [], saveto, summary)
    cursor.rnalib = next(batches.rnalib, None)

while cursor.genome is not None:
    nms([], cursor.genome, saveto, summary)
    cursor.genome = next(batches.genome, None)

assert cursor.rnalib is None and cursor.genome is None

saveto.close()

# Remove temporary files
sorted_by_name.rnalib.unlink()
sorted_by_name.genome.unlink()
del sorted_by_name

# Calculate the total number of reads
total = 0
with pysam.FastxFile(parser.reads, persist=False) as fastq:
    for r in fastq:
        total += 1
summary["Unmapped"] = total - sum(summary.values())
assert summary["Unmapped"] >= 0

# Sort & Index
tmpsort = parser.alignments.with_suffix(f".sorted.bam")
pysam.sort(
    "-o", tmpsort.as_posix(), "-@", str(parser.threads), parser.alignments.as_posix()
)
tmpsort.rename(parser.alignments)

pysam.index("-@", str(parser.threads), parser.alignments.as_posix())

# Write summary
print(f"Writing NMS summary to {parser.logs}....")
sample = parser.logs.name.replace("_nms.tsv", "")

items = sorted(summary.items())
with open(parser.logs, "w") as saveto:
    # Header
    saveto.write("Sample")
    for k, _ in items:
        saveto.write(f"\t{k}")
    saveto.write("\n")

    # Data
    saveto.write(f"{sample}")
    for _, v in items:
        saveto.write(f"\t{v}")
    saveto.write("\n")
