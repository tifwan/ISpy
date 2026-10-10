#!/usr/bin/env python3
"""Combine IS reference databases, add novel IS elements found by ISEScan.

Steps
  1. Concatenate the reference FASTAs (in the order given; earlier files win
     ties, so put ISFinder first). Drop exact-duplicate sequences (either
     strand) and make IDs unique (ISFinder IS.fna repeats some names).
  2. Pool ISEScan *.is.fna sequences and collapse redundant copies with
     cd-hit-est (default 95% identity, 90% coverage of the shorter seq).
  3. BLAST the cluster representatives against the combined database. A
     representative is "known" if any hit has >= --identity % identity and
     covers >= --qcov % of the query; otherwise it is novel. (95% nt identity
     is the ISfinder threshold for calling two elements the same IS.)
  4. Rename novel elements IS_novel_NNNN_<family> and append them to the
     database. Numbering continues from any IS_novel_ (or older novelIS_) IDs
     already present, so a previous output can be passed back in as --db for
     incremental updates.
Without --isescan/--panisa only step 1 runs (used to build the combined db).

Outputs in --outdir (prefix set by --name):
  <name>.fna                 combined + novel database (makeblastdb'd)
  <name>.ref_duplicates.tsv  reference sequences dropped/renamed in step 1
  <name>.isescan_clusters.tsv  every ISEScan seq -> its representative
  <name>.isescan_vs_db.tsv   best DB hit per representative, known/novel call
  <name>.novel.fna           novel elements only, with new names
  <name>.novel_map.tsv       new name -> original ISEScan header
  <name>.panisa_rescue.tsv   per panISa site: element extracted from assembly?
  <ISEScan file>.renamed.fna  elements renamed to DB/novel IDs (in --renamed-dir,
                              default next to each input)

Optional panISa rescue (--panisa): ISEScan misses some elements. For each
panISa site without a DB hit, the 30 bp IS ends (same as local_db_blast.py uses) are located in
the sample assembly (--assembly-pattern) and the stretch between them is added to the
ISEScan pool, so elements behind "No identity" sites also enter the database.

Requires blastn, makeblastdb and cd-hit-est on PATH (or pass --cdhit).
"""

import argparse
import importlib.util
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

COMP = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def write_fasta(fh, name, seq, width=60):
    fh.write(">%s\n" % name)
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def canonical(seq):
    s = seq.upper()
    rc = s.translate(COMP)[::-1]
    return min(s, rc)


def run(cmd):
    print("+ " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def combine_refs(paths, dup_log):
    """Return list of (id, description, seq) with exact dups removed."""
    seen_seq = {}          # canonical seq -> kept id
    used_ids = defaultdict(int)
    records = []
    with open(dup_log, "w") as log:
        log.write("source\toriginal_header\taction\tkept_as\n")
        for path in paths:
            src = os.path.basename(path)
            for header, seq in read_fasta(path):
                sid = header.split()[0]
                key = canonical(seq)
                if key in seen_seq:
                    log.write("%s\t%s\tdropped_identical_seq\t%s\n" % (src, header, seen_seq[key]))
                    continue
                used_ids[sid] += 1
                new_id = sid if used_ids[sid] == 1 else "%s_v%d" % (sid, used_ids[sid])
                if new_id != sid:
                    log.write("%s\t%s\trenamed_duplicate_id\t%s\n" % (src, header, new_id))
                seen_seq[key] = new_id
                records.append((new_id, header[len(sid):].strip(), seq))
    return records, seen_seq


def parse_cdhit_clusters(clstr):
    """Map every member id -> representative id."""
    rep_of, members = {}, []
    rep = None

    def flush():
        for m in members:
            rep_of[m] = rep

    with open(clstr) as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                if members:
                    flush()
                members, rep = [], None
                continue
            m = re.search(r">(\S+?)\.\.\.", line)
            if not m:
                continue
            members.append(m.group(1))
            if line.rstrip().endswith("*"):
                rep = m.group(1)
    if members:
        flush()
    return rep_of


def load_panisa_helpers():
    """Import read_panisa/get_irl_irr from local_db_blast.py (same directory)."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "local_db_blast.py")
    spec = importlib.util.spec_from_file_location("local_db_blast", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rescue_from_panisa(a, tmp, log_path, refs_db):
    """Cut out the element between the two panISa IS ends in the sample assembly.

    Only sites whose ends match nothing in refs_db (ISFinder_local.py default
    filters, i.e. what would be reported as "No identity") are rescued; ends
    that already hit a known IS can sit inside it, giving a bogus span.
    Returns a list of (header, seq). Every site is logged with a status.
    """
    helpers = load_panisa_helpers()
    files = [open(f) for f in a.panisa]
    merged = helpers.read_panisa(files, a.panisa_suffix)
    for f in files:
        f.close()

    by_sample = defaultdict(list)
    for v in merged:
        by_sample[v["sample"]].append(v)

    L = 30
    rescued = []
    with open(log_path, "w") as log:
        log.write("sample\tsite\tstatus\tcontig\tstart\tend\tstrand\tlength\n")
        for sample, sites in by_sample.items():
            asm = a.assembly_pattern.format(sample=sample)
            site_seqs = helpers.get_irl_irr(sites, L)
            raw = helpers.run_blast(site_seqs, refs_db, a.threads)
            matched = {cp for cp, lines in raw.items() if helpers.filter_blast(lines, 90, 0.001, L, 80)}
            ends = {}
            for v in sites:
                both = site_seqs[v["chrompos"]]
                left, right = both[:L], both[L:]
                site = "%s_%s" % (v["chrom"].split("|")[-2] if "|" in v["chrom"] else v["chrom"], v["lpos"])
                if v["chrompos"] in matched:
                    log.write("%s\t%s\talready_matched\t\t\t\t\t\n" % (sample, site))
                    continue
                if left.count("N") > 10 or right.count("N") > 10:
                    log.write("%s\t%s\tone_sided\t\t\t\t\t\n" % (sample, site))
                    continue
                ends[site] = (left, right)
            if not ends:
                continue
            if not os.path.exists(asm):
                for site in ends:
                    log.write("%s\t%s\tno_assembly\t\t\t\t\t\n" % (sample, site))
                continue

            qfa = os.path.join(tmp, "%s_ends.fa" % sample)
            with open(qfa, "w") as fh:
                for site, (left, right) in ends.items():
                    fh.write(">%s|L\n%s\n>%s|R\n%s\n" % (site, left, site, right))
            res = subprocess.run(
                ["blastn", "-task", "blastn-short", "-query", qfa, "-subject", asm, "-evalue", "1",
                 "-outfmt", "6 qseqid sseqid pident length sstart send", "-max_target_seqs", "100"],
                check=True, capture_output=True, text=True).stdout
            hits = defaultdict(list)     # (site, side) -> [(contig, lo, hi, strand)]
            for line in res.splitlines():
                q, s, pident, length, ss, se = line.split("\t")
                if float(pident) < 95 or int(length) < 28:
                    continue
                site, side = q.rsplit("|", 1)
                ss, se = int(ss), int(se)
                hits[(site, side)].append((s, min(ss, se), max(ss, se), "+" if se > ss else "-"))

            contigs = None
            for site in ends:
                lh, rh = hits.get((site, "L"), []), hits.get((site, "R"), [])
                if not lh or not rh:
                    log.write("%s\t%s\tends_not_found\t\t\t\t\t\n" % (sample, site))
                    continue
                pairs = []
                for lc, llo, lhi, lst in lh:
                    for rc, rlo, rhi, rst in rh:
                        if lc != rc or lst != rst:
                            continue
                        lo, hi = min(llo, rlo), max(lhi, rhi)
                        if a.rescue_min <= hi - lo + 1 <= a.rescue_max:
                            pairs.append((hi - lo + 1, lc, lo, hi, lst))
                if not pairs:
                    log.write("%s\t%s\tno_valid_pair\t\t\t\t\t\n" % (sample, site))
                    continue
                length, contig, lo, hi, strand = min(pairs)
                if contigs is None:
                    contigs = {h.split()[0]: s for h, s in read_fasta(asm)}
                seq = contigs[contig][lo - 1:hi]
                if strand == "-":
                    seq = seq.translate(COMP)[::-1]
                rescued.append(("%s_%s panISa_rescue %s:%d-%d(%s)" % (sample, site, contig, lo, hi, strand), seq))
                log.write("%s\t%s\textracted\t%s\t%d\t%d\t%s\t%d\n" % (sample, site, contig, lo, hi, strand, length))
    return rescued


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", nargs="+", required=True, help="Reference IS FASTA(s), highest priority first")
    p.add_argument("--isescan", nargs="*", default=[], help="ISEScan *.is.fna file(s)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--name", default="IS_combined")
    p.add_argument("--identity", type=float, default=95, help="Min %% identity to call known (default 95)")
    p.add_argument("--qcov", type=float, default=90, help="Min %% query coverage to call known (default 90)")
    p.add_argument("--cluster-id", type=float, default=0.95, help="cd-hit-est -c for ISEScan seqs (default 0.95)")
    p.add_argument("--cluster-cov", type=float, default=0.9, help="cd-hit-est -aS (default 0.9)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--cdhit", default=shutil.which("cd-hit-est") or "cd-hit-est")
    p.add_argument("--panisa", nargs="*", default=[], help="panISa *_panISa.txt file(s) for rescue")
    p.add_argument("--panisa-suffix", default="_panISa.txt", help="Stripped from panISa file names to get sample")
    p.add_argument("--assembly-pattern", default="assembly/{sample}/{sample}.fna",
                   help="Sample assembly path, {sample} is filled in (default assembly/{sample}/{sample}.fna)")
    p.add_argument("--renamed-dir", help="Write <ISEScan file>.renamed.fna here (default: next to each input)")
    p.add_argument("--rescue-min", type=int, default=150, help="Min rescued element length (default 150)")
    p.add_argument("--rescue-max", type=int, default=10000, help="Max rescued element length (default 10000)")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    out = lambda suffix: os.path.join(a.outdir, a.name + suffix)
    tmp = tempfile.mkdtemp(prefix="isdb_", dir=a.outdir)

    # 1. combine references
    refs, seen_seq = combine_refs(a.db, out(".ref_duplicates.tsv"))
    ref_fna = os.path.join(tmp, "refs.fna")
    with open(ref_fna, "w") as fh:
        for sid, _, seq in refs:
            write_fasta(fh, sid, seq)
    n_in = sum(1 for f in a.db for _ in read_fasta(f))
    print("references: %d input -> %d after exact-dup removal" % (n_in, len(refs)), file=sys.stderr)

    # 2. pool + cluster ISEScan sequences (IDs made safe for cd-hit/blast)
    iss, orig, src = {}, {}, {}
    for path in a.isescan:
        for header, seq in read_fasta(path):
            qid = "q%05d" % (len(iss) + 1)
            iss[qid] = seq
            orig[qid] = header
            src[qid] = path
    run(["makeblastdb", "-in", ref_fna, "-dbtype", "nucl", "-out", os.path.join(tmp, "refs")])
    if a.panisa:
        rescued = rescue_from_panisa(a, tmp, out(".panisa_rescue.tsv"), os.path.join(tmp, "refs"))
        print("panISa rescue: %d elements extracted from assemblies" % len(rescued), file=sys.stderr)
        for header, seq in rescued:
            qid = "q%05d" % (len(iss) + 1)
            iss[qid] = seq
            orig[qid] = header
    else:
        with open(out(".panisa_rescue.tsv"), "w") as log:
            log.write("sample\tsite\tstatus\tcontig\tstart\tend\tstrand\tlength\n")
    pooled = os.path.join(tmp, "isescan.fna")
    with open(pooled, "w") as fh:
        for qid, seq in iss.items():
            write_fasta(fh, qid, seq)
    clustered = os.path.join(tmp, "isescan.cdhit.fna")
    rep_of = {}
    if iss:                              # cd-hit-est fails on an empty file
        run([a.cdhit, "-i", pooled, "-o", clustered, "-c", str(a.cluster_id), "-aS", str(a.cluster_cov),
             "-n", "10", "-G", "0", "-r", "1", "-d", "0", "-M", "0", "-T", str(a.threads)])
        rep_of = parse_cdhit_clusters(clustered + ".clstr")
    reps = sorted(set(rep_of.values()))
    with open(out(".isescan_clusters.tsv"), "w") as fh:
        fh.write("isescan_header\trepresentative_header\n")
        for qid in iss:
            fh.write("%s\t%s\n" % (orig[qid], orig[rep_of[qid]]))
    print("ISEScan: %d sequences -> %d clusters" % (len(iss), len(reps)), file=sys.stderr)

    # 3. BLAST representatives vs combined refs
    blast_out = os.path.join(tmp, "reps_vs_refs.tsv")
    open(blast_out, "w").close()
    if reps:
        run(["blastn", "-query", clustered, "-db", os.path.join(tmp, "refs"), "-out", blast_out,
             "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore qcovs",
             "-evalue", "1e-10", "-max_target_seqs", "50", "-num_threads", str(a.threads)])
    best = {}
    with open(blast_out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, pident, qcov, bits = f[0], float(f[2]), float(f[10]), float(f[9])
            known = pident >= a.identity and qcov >= a.qcov
            score = (known, bits)
            if q not in best or score > best[q][0]:
                best[q] = (score, f)

    # 4. name novel reps and write outputs
    existing = [int(m.group(1)) for sid, _, _ in refs for m in [re.match(r"(?:IS_novel|novelIS)_(\d+)_", sid)] if m]
    n_next = max(existing, default=0) + 1
    novel = []
    rep_name = {}                        # rep -> (name, call, pident, qcov)
    with open(out(".isescan_vs_db.tsv"), "w") as fh:
        fh.write("isescan_header\tcluster_size\tcall\tnew_name\tbest_hit\tpident\tqcov\tqlen\tslen\tevalue\n")
        sizes = defaultdict(int)
        for r in rep_of.values():
            sizes[r] += 1
        for q in reps:
            key = canonical(iss[q])
            hit = best.get(q)
            f = hit[1] if hit else None
            if key in seen_seq:          # identical to a ref sequence already
                call, new = "known", ""
                rep_name[q] = (seen_seq[key], call, "100", "100")
            elif hit and hit[0][0]:
                call, new = "known", ""
                rep_name[q] = (f[1], call, f[2], f[10])
            else:
                fam = orig[q].split()[1].rsplit("_", 1)[0] if len(orig[q].split()) > 1 else "unknown"
                new = "IS_novel_%04d_%s" % (n_next, fam)
                n_next += 1
                call = "novel"
                novel.append((new, q))
                rep_name[q] = (new, call, f[2] if f else "", f[10] if f else "")
            if f:
                fh.write("\t".join([orig[q], str(sizes[q]), call, new, f[1], f[2], f[10], f[4], f[5], f[8]]) + "\n")
            else:
                fh.write("\t".join([orig[q], str(sizes[q]), call, new, "", "", "", str(len(iss[q])), "", ""]) + "\n")

    with open(out(".novel.fna"), "w") as nf, open(out(".novel_map.tsv"), "w") as mf:
        mf.write("new_name\tisescan_header\tlength\tcluster_size\n")
        for new, q in novel:
            write_fasta(nf, "%s %s" % (new, orig[q]), iss[q])
            mf.write("%s\t%s\t%d\t%d\n" % (new, orig[q], len(iss[q]), sizes[q]))

    with open(out(".fna"), "w") as fh:
        for sid, desc, seq in refs:
            write_fasta(fh, (sid + " " + desc).strip(), seq)
        for new, q in novel:
            write_fasta(fh, "%s %s" % (new, orig[q]), iss[q])
    run(["makeblastdb", "-in", out(".fna"), "-dbtype", "nucl", "-out", out(".fna")])

    # renamed copy of each ISEScan FASTA, in --renamed-dir or next to the original
    if a.renamed_dir:
        os.makedirs(a.renamed_dir, exist_ok=True)
    for path in a.isescan:
        renamed = re.sub(r"\.fna$", "", path) + ".renamed.fna"
        if a.renamed_dir:
            renamed = os.path.join(a.renamed_dir, os.path.basename(renamed))
        with open(renamed, "w") as fh:
            for qid in (q for q in iss if src.get(q) == path):
                name, call, pident, qcov = rep_name[rep_of[qid]]
                write_fasta(fh, "%s %s call=%s pident=%s qcov=%s" % (name, orig[qid], call, pident, qcov), iss[qid])
        print("renamed: " + renamed, file=sys.stderr)

    shutil.rmtree(tmp)
    print("novel: %d / %d clusters; database: %s (%d seqs)"
          % (len(novel), len(reps), out(".fna"), len(refs) + len(novel)), file=sys.stderr)


if __name__ == "__main__":
    main()
