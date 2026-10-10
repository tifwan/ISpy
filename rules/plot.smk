# Step 8: Plot insertion sites and their neighbouring genes

rule plot_IS_sites:
    input:
        sites = f"{RES}/{{sample}}/annotation/{{sample}}_{{db}}_sites.bed",
        genes = f"{RES}/reference/genes.bed",
        fai = REF + ".fai",
        families = config["isosdb_family_annotation"]
    output:
        neighbourhood = f"{RES}/{{sample}}/plots/{{sample}}_{{db}}_neighbourhood.pdf",
        overview = f"{RES}/{{sample}}/plots/{{sample}}_{{db}}_overview.pdf",
        compressed = f"{RES}/{{sample}}/plots/{{sample}}_{{db}}_compressed.pdf"
    params:
        flank = config["plot"]["flank"]
    log:
        f"{LOGS}/{{sample}}/plots/{{sample}}_{{db}}_plot.log"
    conda:
        "../envs/r_plot.yaml"
    script:
        "../scripts/plot_IS_sites.R"


rule compare_dbs:
    input:
        isfinder = f"{RES}/{{sample}}/blast/{{sample}}_ISFinder.tsv",
        isosdb = f"{RES}/{{sample}}/blast/{{sample}}_ISOSDB.tsv",
        families = config["isosdb_family_annotation"]
    output:
        table = f"{RES}/{{sample}}/compare/{{sample}}_ISFinder_vs_ISOSDB.tsv",
        venn = f"{RES}/{{sample}}/plots/{{sample}}_ISFinder_vs_ISOSDB_venn.pdf"
    log:
        f"{LOGS}/{{sample}}/plots/{{sample}}_compare_dbs.log"
    conda:
        "../envs/r_plot.yaml"
    script:
        "../scripts/compare_dbs.R"
