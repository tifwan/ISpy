# ISpy, workflow 1 of 3: IS detection on the existing databases
#   trim -> align -> dedup -> panISa -> local BLAST (ISFinder / ISOSDB / combined)
#   -> genes near each insertion (bedtools) -> plots; plus ISEScan on the assemblies
# Based on ISScreener (Kyle Gontjes): https://github.com/kylegontjes/ISScreener
# Run order: ISpy.smk -> update_db.smk -> reblast.smk

configfile: "config/config.yaml"

include: "rules/common.smk"


rule all:
    input:
        expand(f"{RES}/{{sample}}/align_reads/{{sample}}_final.bam", sample=SAMPLES),
        expand(f"{RES}/{{sample}}/panISa/{{sample}}_panISa.txt", sample=SAMPLES),
        expand(f"{RES}/isescan/{{sample}}/{{sample}}.is.fna", sample=SAMPLES),
        db_targets(config["blast_dbs"]),
        compare_targets()


include: "rules/preprocess.smk"
include: "rules/panisa.smk"
include: "rules/blast.smk"
include: "rules/isescan.smk"
include: "rules/annotate.smk"
include: "rules/plot.smk"
