#!/usr/bin/env python3
"""bedtools window/closest output -> one row per (insertion site, gene).

Inputs are the sites BED (make_sites_bed.py, 9 columns) and the bedtools output
of sites vs genes BED (gff_to_bed.py, 9 columns); `bedtools closest -d` adds a
distance column, which is ignored and recomputed here.

distance: bp between site and gene (0 = the IS sits inside the gene)
position: inside_gene, or upstream/downstream relative to the gene's strand
Sites with no gene reported by bedtools get one row with position "no_gene".

  python summarize_sites.py --sites s.bed --hits s_bedtools.tsv -o s_IS_genes.tsv
"""
import argparse
from collections import defaultdict

N = 9  # columns in each BED

HEADER = ["sample", "site", "contig", "start", "end", "Potential_IS", "Alignment",
          "locus_tag", "gene", "product", "feature_type", "gene_start", "gene_end",
          "gene_strand", "distance", "position"]


def relate(s0, e, g0, g1, strand):
    """Site [s0, e) vs gene [g0, g1), both BED half-open."""
    if s0 < g1 and g0 < e:
        return 0, "inside_gene"
    left_of_gene = e <= g0
    gap = g0 - e + 1 if left_of_gene else s0 - g1 + 1
    upstream = left_of_gene if strand != "-" else not left_of_gene
    return gap, "upstream" if upstream else "downstream"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sites", required=True)
    p.add_argument("--hits", required=True)
    p.add_argument("-o", "--output", required=True)
    a = p.parse_args()

    hits = defaultdict(list)
    with open(a.hits) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 * N:
                continue
            gene = f[N:2 * N]
            if gene[0] == "." or gene[1] == "-1":       # closest: no gene on this contig
                continue
            hits[tuple(f[:N])].append(gene)

    with open(a.sites) as fh, open(a.output, "w") as out:
        out.write("\t".join(HEADER) + "\n")
        for line in fh:
            site = line.rstrip("\n").split("\t")
            if len(site) < N:
                continue
            contig, s0, e, name, _, _, sample, is_name, aln = site[:N]
            s0, e = int(s0), int(e)
            base = [sample, name, contig, str(s0 + 1), str(e), is_name, aln]
            rows = []
            for g in hits.get(tuple(site[:N]), []):
                g0, g1 = int(g[1]), int(g[2])
                dist, pos = relate(s0, e, g0, g1, g[5])
                rows.append((dist, g0, [g[3], g[6], g[7], g[8], str(g0 + 1), str(g1), g[5], str(dist), pos]))
            if not rows:
                rows.append((0, 0, ["."] * 8 + ["no_gene"]))
            for _, _, r in sorted(rows, key=lambda x: (x[0], x[1])):
                out.write("\t".join(base + r) + "\n")


if __name__ == "__main__":
    main()
