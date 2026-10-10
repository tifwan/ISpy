#!/usr/bin/env python3
"""
Local replacement for panISa's ISFinder_search.py.

Same logic, same defaults, same output table -- but BLASTs against a local
ISfinder BLAST database instead of the www-is.biotoul.fr web server.

  makeblastdb -in IS.fna -dbtype nucl -out IS_db
  python ISFinder_local.py -d IS_db -o sample_ISFinder.txt sample_panISa.txt
"""
import argparse, subprocess, sys, tempfile, os

HEADER = ["chrom", "lpos", "lclipped", "dr", "rpos", "rclipped", "ir", "irl", "irr"]

def read_panisa(files, remove):
    merged = []
    for fi in files:
        fi.readline()                       # discard header
        for line in fi:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            tab = line.split("\t")
            if len(tab) != 9:
                raise Exception("Not a panISa output or incorrect\n" + line)
            r = dict(zip(HEADER, tab))
            name = os.path.basename(fi.name)
            if name.endswith(remove):
                name = name[: -len(remove)]
            r["sample"] = name
            r["chrompos"] = r["chrom"] + "_" + r["lpos"]
            merged.append(r)
    return merged

def get_irl_irr(merged, length):
    seqs = {}
    for v in merged:
        irl = v["irl"][:length]
        irr = v["irr"][-length:] if v["irr"] else ""
        irr = irr + "N" * (length - len(irr))
        irl = "N" * (length - len(irl)) + irl
        cp, cur = v["chrompos"], irl + irr
        if cp not in seqs or cur.count("N") < seqs[cp].count("N"):
            seqs[cp] = cur
    return seqs

def run_blast(seqs, db, threads):
    """One blastn call for all queries; same params as the ISFinder web form."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for cp, seq in seqs.items():
            fh.write(">%s\n%s\n" % (cp, seq))
        qpath = fh.name
    cmd = ["blastn", "-query", qpath, "-db", db,
           "-word_size", "11", "-gapopen", "5", "-gapextend", "2",
           "-evalue", "10.0", "-outfmt", "6", "-num_threads", str(threads),
           "-max_target_seqs", "50"]
    out = subprocess.run(cmd, capture_output=True, text=True)
    os.unlink(qpath)
    if out.returncode != 0:
        sys.exit("blastn failed:\n" + out.stderr)
    per_query = {cp: [] for cp in seqs}
    for line in out.stdout.rstrip("\n").split("\n"):
        if not line:
            continue
        per_query.setdefault(line.split("\t")[0], []).append(line)
    return per_query

def filter_blast(lines, identity, evalue, length, alignment):
    cols = ["Query","Subject_ID","%identity","alignment_length","mismatches",
            "gap_opens","q.start","q.end","s.start","s.end","evalue","bit_score"]
    keep = []
    for line in lines:
        b = line.split("\t")
        if len(b) != 12:
            sys.stderr.write(line + "\n"); continue
        r = dict(zip(cols, b))
        aln = int(r["q.end"]) - int(r["q.start"])
        if (float(r["evalue"]) <= evalue and float(r["%identity"]) >= identity
                and aln >= (alignment * length) / 100):
            keep.append(line)
    keep.sort(key=lambda l: float(l.split("\t")[11]), reverse=True)
    return keep

def write_result(output, merged, seqs, length, results):
    output.write("Sample\tChromosome\tStart_Position\tStop_Position\t"
                 "Potential_sequence\tPotential_IS\tAlignment\n")
    for v in merged:
        cp = v["chrompos"]
        seq = seqs[cp]
        row = [v["sample"], v["chrom"], v["rpos"], v["lpos"],
               seq[:length] + "-" + seq[-length:]]
        hits = results.get(cp, [])
        if len(hits) == 0:
            row += ["No identity", ""]
        elif len(hits) == 1:
            row += [hits[0].split("\t")[1], "one side"]
        else:
            first = hits[0].split("\t")[1]
            row += [first, "two side" if first == hits[1].split("\t")[1] else "one side"]
        output.write("\t".join(row) + "\n")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="local ISFinder search from panISa output")
    p.add_argument("-d", "--db", required=True, help="local BLAST db prefix (e.g. IS_db)")
    p.add_argument("-o", "--output", nargs="?", type=argparse.FileType("wt"), default=sys.stdout)
    p.add_argument("-r", "--remove", default=".txt")
    p.add_argument("-l", "--length", type=int, default=30)
    p.add_argument("-i", "--identity", type=int, default=90)
    p.add_argument("-e", "--evalue", type=float, default=0.001)
    p.add_argument("-a", "--alignment", type=float, default=80)
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("file", nargs="+", type=argparse.FileType("r"))
    a = p.parse_args()

    merged = read_panisa(a.file, a.remove)
    seqs = get_irl_irr(merged, a.length)
    raw = run_blast(seqs, a.db, a.threads)
    results = {cp: filter_blast(v, a.identity, a.evalue, a.length, a.alignment)
               for cp, v in raw.items()}
    write_result(a.output, merged, seqs, a.length, results)
