#!/usr/bin/env python
import argparse
import pysam
from pathlib import Path
from math import log10

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


parser = argparse.ArgumentParser(
    description="Fix alignment flags & HN tag after UMI deduplication",
)
parser.add_argument("threads", type=int, help="Number of threads to use")
parser.add_argument("bam", type=Path, help="Input bam file")
parser.add_argument("saveto", type=Path, help="Where to save results")

# args = ["64", "reads.fastq.gz", "smRNA.bam", "genome.bam", "stats.tsv", "nms.bam"]
# parser = parser.parse_args(args)

parser = parser.parse_args()

# Sort by name input files
print("Sorting input reads by name....")
bam = parser.bam.with_suffix(".sorted-by-name.bam")
pysam.sort("-@", str(parser.threads), "-o", bam.as_posix(), "-n", parser.bam.as_posix())


with pysam.AlignmentFile(bam, 'rb') as template:
    saveto = pysam.AlignmentFile(parser.saveto, 'wb', header=template.header)

batches = batch_reads(bam)

for batch in batches:
    # Update tags/flags
    mapq = 255 if len(batch) == 1 else int(-10 * log10(1 - 1 / len(batch)))
    for read in batch:
        read.set_tag("NH", len(batch))
        read.mapping_quality = mapq
        read.is_secondary = True
    batch[0].is_secondary = False

    # Write back the result
    for read in batch:
        saveto.write(read)

saveto.close()
bam.unlink()

# Sort & Index
tmpsort = parser.saveto.with_suffix(f".sorted.bam")
pysam.sort(
    "-o", tmpsort.as_posix(), "-@", str(parser.threads), parser.saveto.as_posix()
)
tmpsort.rename(parser.saveto)

pysam.index("-@", str(parser.threads), parser.saveto.as_posix())
