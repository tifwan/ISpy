# Step 7: Genes at / near each insertion site (bedtools)
#   method: window  -> every gene within `distance` bp of the site
#   method: closest -> the `n_closest` nearest genes, any distance

ANNOT = config["annotation"]

if ANNOT["method"] == "window":
    BEDTOOLS_CMD = f"bedtools window -a {{input.sites}} -b {{input.genes}} -w {ANNOT['distance']}"
elif ANNOT["method"] == "closest":
    BEDTOOLS_CMD = f"bedtools closest -a {{input.sites}} -b {{input.genes}} -k {ANNOT['n_closest']} -d -t all"
else:
    raise ValueError("annotation: method must be 'window' or 'closest'")


rule ref_genes_bed:
    input:
        gff = config["reference_gff"]
    output:
        f"{RES}/reference/genes.bed"
    params:
        types = " ".join(ANNOT["feature_types"])
    log:
        f"{LOGS}/reference/ref_genes_bed.log"
    conda:
        "../envs/blast.yaml"
    shell:
        '(python "{SCRIPTS}/gff_to_bed.py" {input.gff} --types {params.types} '
        "| LC_ALL=C sort -k1,1 -k2,2n > {output}) 2> {log}"


rule sites_bed:
    input:
        f"{RES}/{{sample}}/blast/{{sample}}_{{db}}.tsv"
    output:
        f"{RES}/{{sample}}/annotation/{{sample}}_{{db}}_sites.bed"
    log:
        f"{LOGS}/{{sample}}/annotation/{{sample}}_{{db}}_sites_bed.log"
    conda:
        "../envs/blast.yaml"
    shell:
        '(python "{SCRIPTS}/make_sites_bed.py" {input} '
        "| LC_ALL=C sort -k1,1 -k2,2n > {output}) 2> {log}"


rule nearby_genes:
    input:
        sites = f"{RES}/{{sample}}/annotation/{{sample}}_{{db}}_sites.bed",
        genes = f"{RES}/reference/genes.bed"
    output:
        f"{RES}/{{sample}}/annotation/{{sample}}_{{db}}_bedtools_{ANNOT['method']}.tsv"
    log:
        f"{LOGS}/{{sample}}/annotation/{{sample}}_{{db}}_bedtools.log"
    conda:
        "../envs/blast.yaml"
    shell:
        BEDTOOLS_CMD + " > {output} 2> {log}"


rule summarize_sites:
    input:
        sites = f"{RES}/{{sample}}/annotation/{{sample}}_{{db}}_sites.bed",
        hits = f"{RES}/{{sample}}/annotation/{{sample}}_{{db}}_bedtools_{ANNOT['method']}.tsv"
    output:
        f"{RES}/{{sample}}/annotation/{{sample}}_{{db}}_IS_genes.tsv"
    log:
        f"{LOGS}/{{sample}}/annotation/{{sample}}_{{db}}_summarize.log"
    conda:
        "../envs/blast.yaml"
    shell:
        'python "{SCRIPTS}/summarize_sites.py" --sites {input.sites} --hits {input.hits} '
        "-o {output} 2> {log}"
