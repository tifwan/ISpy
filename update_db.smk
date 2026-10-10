# ISpy, workflow 2 of 3: build an updated IS database
#   base databases (e.g. ISFinder + ISOSDB) + ISEScan elements from every sample
#   -> deduplicate -> elements matching nothing become IS_novel_* -> new database
# Needs the panISa results from ISpy.smk when update_db: panisa_rescue is true.

configfile: "config/config.yaml"

include: "rules/common.smk"

RESCUE = UPDATE.get("panisa_rescue", True)
OUT = f"{UPDATE['outdir']}/{NEW_DB}"


rule all:
    input:
        f"{OUT}.fna"


include: "rules/blast.smk"
include: "rules/isescan.smk"


rule update_db:
    input:
        base = [DB_FASTA[d] for d in UPDATE["base"]],
        isescan = expand(f"{RES}/isescan/{{sample}}/{{sample}}.is.fna", sample=SAMPLES),
        panisa = expand(f"{RES}/{{sample}}/panISa/{{sample}}_panISa.txt", sample=SAMPLES) if RESCUE else [],
        assemblies = expand(config["assembly"], sample=SAMPLES) if RESCUE else []
    output:
        fna = f"{OUT}.fna",
        novel = f"{OUT}.novel.fna",
        novel_map = f"{OUT}.novel_map.tsv",
        ref_duplicates = f"{OUT}.ref_duplicates.tsv",
        clusters = f"{OUT}.isescan_clusters.tsv",
        vs_db = f"{OUT}.isescan_vs_db.tsv",
        rescue = f"{OUT}.panisa_rescue.tsv",
        renamed = directory(f"{UPDATE['outdir']}/renamed_isescan")
    params:
        outdir = UPDATE["outdir"],
        name = NEW_DB,
        panisa = lambda wildcards, input: ("--panisa " + " ".join(input.panisa)) if RESCUE else "",
        assembly = config["assembly"],
        identity = UPDATE["identity"],
        qcov = UPDATE["qcov"],
        cluster_id = UPDATE["cluster_id"],
        cluster_cov = UPDATE["cluster_cov"]
    threads: THREADS
    log:
        f"{LOGS}/update_db/{NEW_DB}.log"
    conda:
        "envs/blast.yaml"
    shell:
        'python "{SCRIPTS}/update_IS_db.py" --db {input.base} --isescan {input.isescan} {params.panisa} '
        "--assembly-pattern '{params.assembly}' --outdir {params.outdir} --name {params.name} "
        "--renamed-dir {output.renamed} --identity {params.identity} --qcov {params.qcov} "
        "--cluster-id {params.cluster_id} --cluster-cov {params.cluster_cov} "
        "--threads {threads} &> {log}"
