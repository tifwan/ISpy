#!/usr/bin/env python3
"""Reference GFF3 (Prokka or NCBI) -> BED of genes for bedtools.

Columns: contig  start(0-based)  end  locus_tag  .  strand  gene  product  feature_type
Reading stops at ##FASTA. Contig names such as gi|661922017|gb|CP008827.1| are
shortened to the accession (CP008827.1), as in make_sites_bed.py.
Output is unsorted; pipe through `LC_ALL=C sort -k1,1 -k2,2n`.

  python gff_to_bed.py KPNIH1.gff --types CDS tRNA rRNA > genes.bed
"""
import argparse
import sys
from urllib.parse import unquote


def norm_seqid(seqid):
    parts = [p for p in seqid.split("|") if p]
    return parts[-1] if parts else seqid


def parse_attrs(col):
    attrs = {}
    for kv in col.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            attrs[k] = unquote(v).replace("\t", " ")
    return attrs


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("gff")
    p.add_argument("--types", nargs="+", default=["CDS", "tRNA", "rRNA"],
                   help="GFF feature types to keep (default CDS tRNA rRNA)")
    a = p.parse_args()
    types = set(a.types)

    n = 0
    with open(a.gff) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] not in types:
                continue
            at = parse_attrs(f[8])
            locus = at.get("locus_tag") or at.get("ID") or "."
            print("\t".join([norm_seqid(f[0]), str(int(f[3]) - 1), f[4], locus, ".", f[6],
                             at.get("gene", "."), at.get("product", "."), f[2]]))
            n += 1
    sys.stderr.write("%d features of type %s\n" % (n, ", ".join(sorted(types))))
    if n == 0:
        sys.exit("No features of type %s in %s" % (", ".join(sorted(types)), a.gff))


if __name__ == "__main__":
    main()
