# Step 6: BLAST panISa IS ends against a local IS database
# (replaces ISScreener's ISFinder_search.py, which needs the ISfinder web server)

LB = config["local_blast"]


if config["databases"].get("combined") == "auto":
    # ISFinder + ISOSDB with exact duplicates removed (no ISEScan input)
    rule build_combined_db:
        input:
            [DB_FASTA[d] for d in config["combined_from"]]
        output:
            fna = COMBINED_AUTO,
            duplicates = f"{RES}/combined_db/IS_combined.ref_duplicates.tsv"
        params:
            outdir = os.path.dirname(COMBINED_AUTO)
        log:
            f"{LOGS}/combined_db/build_combined_db.log"
        conda:
            "../envs/blast.yaml"
        shell:
            'python "{SCRIPTS}/update_IS_db.py" --db {input} '
            "--outdir {params.outdir} --name IS_combined &> {log}"


rule make_blastdb:
    input:
        fna = lambda wildcards: DB_FASTA[wildcards.db]
    output:
        f"{RES}/blastdb/{{db}}/{{db}}.nin"
    params:
        prefix = lambda wildcards: blastdb_prefix(wildcards.db)
    log:
        f"{LOGS}/blastdb/{{db}}_makeblastdb.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "makeblastdb -in {input.fna} -dbtype nucl -out {params.prefix} &> {log}"


rule local_db_blast:
    input:
        panISa_out = f"{RES}/{{sample}}/panISa/{{sample}}_panISa.txt",
        db = f"{RES}/blastdb/{{db}}/{{db}}.nin"
    output:
        f"{RES}/{{sample}}/blast/{{sample}}_{{db}}.tsv"
    params:
        prefix = lambda wildcards: blastdb_prefix(wildcards.db),
        length = LB["length"],
        identity = LB["identity"],
        evalue = LB["evalue"],
        alignment = LB["alignment"]
    threads: THREADS
    log:
        f"{LOGS}/{{sample}}/blast/{{sample}}_{{db}}.log"
    conda:
        "../envs/blast.yaml"
    shell:
        'python "{SCRIPTS}/local_db_blast.py" -d {params.prefix} -r _panISa.txt '
        "-l {params.length} -i {params.identity} -e {params.evalue} -a {params.alignment} "
        "-t {threads} -o {output} {input.panISa_out} 2> {log}"


# One table per database across all samples (like create_master_ISFinder_report.sh)
rule master_tables:
    input:
        blast = expand(f"{RES}/{{sample}}/blast/{{sample}}_{{db}}.tsv", sample=SAMPLES, allow_missing=True),
        genes = expand(f"{RES}/{{sample}}/annotation/{{sample}}_{{db}}_IS_genes.tsv", sample=SAMPLES, allow_missing=True)
    output:
        blast = f"{RES}/{PREFIX}_{{db}}_master.tsv",
        genes = f"{RES}/{PREFIX}_{{db}}_IS_genes_master.tsv"
    shell:
        """
        awk 'FNR==1 && NR!=1 {{next}} 1' {input.blast} > {output.blast}
        awk 'FNR==1 && NR!=1 {{next}} 1' {input.genes} > {output.genes}
        """
