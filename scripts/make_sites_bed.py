#!/usr/bin/env python3
"""local_db_blast.py table -> BED of IS insertion sites.

Columns: contig  start(0-based)  end  site  .  .  sample  Potential_IS  Alignment
The site spans panISa's Start/Stop positions (the direct repeat). Contig names
such as gi|661922017|gb|CP008827.1| are shortened to the accession, as in
gff_to_bed.py. Output is unsorted; pipe through `LC_ALL=C sort -k1,1 -k2,2n`.

  python make_sites_bed.py sample_ISFinder.tsv > sample_sites.bed
"""
import csv
import sys


def norm_seqid(seqid):
    parts = [p for p in seqid.split("|") if p]
    return parts[-1] if parts else seqid


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: make_sites_bed.py <local_db_blast table>")
    with open(sys.argv[1]) as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            a, b = int(r["Start_Position"]), int(r["Stop_Position"])
            start, end = min(a, b), max(a, b)
            contig = norm_seqid(r["Chromosome"])
            print("\t".join([contig, str(start - 1), str(end), "%s_%d" % (contig, end), ".", ".",
                             r["Sample"], r["Potential_IS"] or ".", r["Alignment"] or "."]))


if __name__ == "__main__":
    main()
