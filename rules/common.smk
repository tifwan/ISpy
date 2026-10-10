# Shared setup for ISpy.smk, update_db.smk and reblast.smk.
# Defines paths, samples and the IS database lookup; contains no rules.

import os
import re
import sys
import pandas as pd

PREFIX = config["prefix"]
RES = f"results/{PREFIX}"
LOGS = f"logs/{PREFIX}"
# relative path, so shell commands work even when the project path contains spaces
SCRIPTS = os.path.relpath(os.path.join(workflow.basedir, "scripts"))

SAMPLES = list(pd.read_csv(config["samples"], sep="\t")["sample_id"].astype(str))
REF = config["reference_genome"]
THREADS = config.get("threads", 8)

# ---- IS databases ---------------------------------------------------------
# name -> FASTA. "auto" for combined = built from config["combined_from"]
# by update_IS_db.py (exact-duplicate removal). The updated database from
# update_db.smk is always available under config["update_db"]["name"].
UPDATE = config["update_db"]
NEW_DB = UPDATE["name"]
COMBINED_AUTO = f"{RES}/combined_db/IS_combined.fna"

DB_FASTA = {name: (COMBINED_AUTO if path == "auto" else path)
            for name, path in config["databases"].items()}
DB_FASTA[NEW_DB] = os.path.join(UPDATE["outdir"], NEW_DB + ".fna")

for db in config["blast_dbs"] + UPDATE["base"]:
    if db not in DB_FASTA:
        raise ValueError(f"database '{db}' is not defined under 'databases' in the config")

# A FASTA that stops mid-sequence (no final newline) is usually an incomplete copy
for name, path in DB_FASTA.items():
    if os.path.isfile(path) and os.path.getsize(path) > 0:
        with open(path, "rb") as fh:
            fh.seek(-1, os.SEEK_END)
            if fh.read(1) != b"\n":
                sys.stderr.write(f"WARNING: {path} ({name}) does not end with a newline; "
                                 "the file may be truncated.\n")


def blastdb_prefix(db):
    return f"{RES}/blastdb/{db}/{db}"


wildcard_constraints:
    sample = "|".join(re.escape(s) for s in SAMPLES),
    db = "|".join(re.escape(d) for d in DB_FASTA)


# ---- targets --------------------------------------------------------------
def db_targets(dbs):
    """BLAST tables, gene annotation, plots and master tables for each db."""
    per_sample = [
        f"{RES}/{{sample}}/blast/{{sample}}_{{db}}.tsv",
        f"{RES}/{{sample}}/annotation/{{sample}}_{{db}}_IS_genes.tsv",
        f"{RES}/{{sample}}/plots/{{sample}}_{{db}}_neighbourhood.pdf",
        f"{RES}/{{sample}}/plots/{{sample}}_{{db}}_overview.pdf",
        f"{RES}/{{sample}}/plots/{{sample}}_{{db}}_compressed.pdf",
    ]
    targets = expand(per_sample, sample=SAMPLES, db=dbs)
    targets += expand(f"{RES}/{PREFIX}_{{db}}_master.tsv", db=dbs)
    targets += expand(f"{RES}/{PREFIX}_{{db}}_IS_genes_master.tsv", db=dbs)
    return targets


def compare_targets():
    """ISFinder vs ISOSDB family comparison, only when both were searched."""
    if {"ISFinder", "ISOSDB"} <= set(config["blast_dbs"]):
        return expand([f"{RES}/{{sample}}/compare/{{sample}}_ISFinder_vs_ISOSDB.tsv",
                       f"{RES}/{{sample}}/plots/{{sample}}_ISFinder_vs_ISOSDB_venn.pdf"],
                      sample=SAMPLES)
    return []
