## Download IS.fna https://github.com/thanhleviet/ISfinder-sequences/blob/master/IS.fna

## update ISFinder_Search.py 
### updated code  - saved the script below as ISFinder_search_local.py 


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## Local-BLAST replacement for ISFinder_search.py (panISa, Valot et al., GPL).
## The original posted queries to the ISfinder web BLAST, which is no longer
## reachable. This version BLASTs the same IRL/IRR query against a local copy
## of the ISfinder sequences (e.g. IS.fna from thanhleviet/ISfinder-sequences)
## and writes the same 7-column output.
##
## Requires BLAST+ (blastn, makeblastdb) on PATH:
##   conda install -c conda-forge -c bioconda blast

"""Search panISa output for IS homology against a local ISfinder database"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

HEADER = ["chrom", "lpos", "lclipped", "dr", "rpos", "rclipped", "ir", "irl", "irr"]


def build_parser():
    desc = "Search IS homology in a local ISfinder database from panISa output"
    p = argparse.ArgumentParser(prog="ISFinder_search_local.py", description=desc,
                                usage="%(prog)s [options] file [file ...]")
    p.add_argument("-d", "--database", default="IS.fna",
                   help="ISfinder sequences in FASTA format (default: IS.fna)")
    p.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout,
                   help="Output table (default: stdout)")
    p.add_argument("-r", "--remove", default=".txt",
                   help="Suffix removed from file names to get sample names (default: .txt)")
    p.add_argument("-l", "--length", type=int, default=30,
                   help="Length of IRL and IRR taken from each end (default: 30)")
    p.add_argument("-i", "--identity", type=float, default=90,
                   help="Minimum percent identity (default: 90)")
    p.add_argument("-e", "--evalue", type=float, default=0.001,
                   help="Maximum e-value (default: 1E-3)")
    p.add_argument("-a", "--alignment", type=float, default=80,
                   help="Minimum alignment length as percent of --length (default: 80)")
    p.add_argument("-t", "--threads", type=int, default=1,
                   help="BLAST threads (default: 1)")
    p.add_argument("file", nargs="+", type=argparse.FileType("r"),
                   help="panISa result file(s)")
    return p


def remove_suffix(name, suffix):
    # The original used str.rstrip(), which strips characters, not a suffix.
    if suffix and name.endswith(suffix):
        return name[:-len(suffix)]
    return name


def read_panisa(files, remove):
    merged = []
    for fi in files:
        fi.readline()  # header
        sample = remove_suffix(os.path.basename(fi.name), remove)
        for line in fi:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            tab = line.split("\t")
            if len(tab) != 9:
                raise ValueError("Not a panISa output or incorrect line:\n" + line)
            rec = dict(zip(HEADER, tab))
            rec["sample"] = sample
            rec["chrompos"] = rec["chrom"] + "_" + rec["lpos"]
            merged.append(rec)
    return merged


def get_irl_irr(merged, length):
    """IRL = first `length` bp of the Left sequence (IS left end),
    IRR = last `length` bp of the Right sequence (IS right end).
    Keep the least-N version per position across samples."""
    seqs = {}
    for v in merged:
        irl = v["irl"][:length].rjust(length, "N")
        irr = v["irr"][-length:].ljust(length, "N")
        cur = irl + irr
        key = v["chrompos"]
        if key not in seqs or cur.count("N") < seqs[key].count("N"):
            seqs[key] = cur
    return seqs


def ensure_db(fasta):
    for tool in ("blastn", "makeblastdb"):
        if shutil.which(tool) is None:
            sys.exit(f"Error: {tool} not found. Install BLAST+: "
                     "conda install -c conda-forge -c bioconda blast")
    if not os.path.exists(fasta):
        sys.exit(f"Error: database FASTA not found: {fasta}")
    built = [fasta + ext for ext in (".nin", ".nal") if os.path.exists(fasta + ext)]
    if not built or os.path.getmtime(built[0]) < os.path.getmtime(fasta):
        sys.stderr.write(f"Building BLAST database from {fasta}...\n")
        r = subprocess.run(["makeblastdb", "-in", fasta, "-dbtype", "nucl"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit("makeblastdb failed:\n" + r.stderr)


def run_blast(seqs, db, evalue, threads):
    """BLAST every query in one call. Returns {chrompos: [hit fields, ...]}."""
    ids = {}
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as q:
        for n, (key, seq) in enumerate(seqs.items()):
            qid = f"q{n}"
            ids[qid] = key
            q.write(f">{qid}\n{seq}\n")
        qpath = q.name
    # Settings mirror the original web search: blastn, word size 11, gap cost 5/2.
    cmd = ["blastn", "-task", "blastn", "-query", qpath, "-db", db,
           "-outfmt", "6", "-evalue", str(evalue), "-word_size", "11",
           "-gapopen", "5", "-gapextend", "2", "-dust", "no",
           "-max_target_seqs", "500", "-num_threads", str(threads)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
    finally:
        os.unlink(qpath)
    if r.returncode != 0:
        sys.exit("blastn failed:\n" + r.stderr)
    hits = {key: [] for key in seqs}
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) == 12:
            hits[ids[f[0]]].append(f)
    return hits


def filter_hits(hits, identity, evalue, length, alignment):
    keep = []
    for f in hits:
        aln = int(f[7]) - int(f[6]) + 1
        if (float(f[10]) <= evalue and float(f[2]) >= identity
                and aln >= alignment * length / 100):
            keep.append(f)
    return keep


def classify(hits, length, alignment):
    """Pick the best IS. Prefer one hit on both ends (IRL in query 1..length,
    IRR in length+1..2*length); otherwise the best single-end hit.
    A hit counts for a side only if it covers enough of that half, so an
    alignment that merely extends a few bases past the midpoint doesn't count."""
    if not hits:
        return "No identity", ""
    min_cov = alignment * length / 100
    subjects = {}
    for f in hits:
        qs, qe = int(f[6]), int(f[7])
        s = subjects.setdefault(f[1], {"sides": set(), "score": 0.0})
        if min(qe, length) - qs + 1 >= min_cov:
            s["sides"].add("left")
        if qe - max(qs, length + 1) + 1 >= min_cov:
            s["sides"].add("right")
        s["score"] += float(f[11])
    both = {k: v for k, v in subjects.items() if {"left", "right"} <= v["sides"]}
    pool = both or subjects
    best = max(pool, key=lambda k: pool[k]["score"])
    return best, ("two side" if both else "one side")


def write_result(output, merged, length, seqs, results):
    output.write("Sample\tChromosome\tStart_Position\tStop_Position\t"
                 "Potential_sequence\tPotential_IS\tAlignment\n")
    for v in merged:
        seq = seqs[v["chrompos"]]
        is_name, side = results[v["chrompos"]]
        output.write("\t".join([v["sample"], v["chrom"], v["rpos"], v["lpos"],
                                seq[:length] + "-" + seq[-length:],
                                is_name, side]) + "\n")


def main():
    args = build_parser().parse_args()
    merged = read_panisa(args.file, args.remove)
    if not merged:
        sys.stderr.write("No insertions in panISa input; nothing to search.\n")
    seqs = get_irl_irr(merged, args.length)
    results = {}
    if seqs:
        ensure_db(args.database)
        hits = run_blast(seqs, args.database, args.evalue, args.threads)
        for key, h in hits.items():
            kept = filter_hits(h, args.identity, args.evalue, args.length, args.alignment)
            results[key] = classify(kept, args.length, args.alignment)
    write_result(args.output, merged, args.length, seqs, results)


if __name__ == "__main__":
    main()


## to run panISa - download panISa.py 

## SCRIPTS is the directory that the panISa.py script and the ISFinder_search_local.py script are located in.
SCRIPTS="/Users/hofskatr/Library/CloudStorage/GoogleDrive-hofstetterkatrina@gmail.com/My Drive/ASM-Hackathon"
S=PCMP_H106     # sample name
R1=PCMP_H106_R1.fastq.gz
R2=PCMP_H106_R2.fastq.gz
REF=KPNIH1.fasta  # the genome you're comparing against

#1. Trim the reads
conda activate trimgalore
trim_galore --paired --cores 4 -o trimmed "$R1" "$R2"

2. Map to the reference with BWA-MEM
bwa index "$REF"                               # once per reference
bwa mem -t 8 "$REF" trimmed/*_val_1.fq.gz trimmed/*_val_2.fq.gz \
  | samtools sort -@ 4 --write-index -o "$S.bam" -

3. Call insertions 
panISa.py "$S.bam" -o "${S}_panisa.tsv"

4. Identify IS elements
python "$SCRIPTS/ISFinder_search_local.py" -d ../IS.fna -r _panisa.tsv \
  -o "${S}_ISfinder.tsv" "${S}_panisa.tsv"

