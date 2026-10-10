# ISpy, workflow 3 of 3: re-BLAST every panISa site against the updated database
# (update_db.smk output), then re-annotate and re-plot.
# Needs the panISa results from ISpy.smk and the database from update_db.smk.

configfile: "config/config.yaml"

include: "rules/common.smk"


rule all:
    input:
        db_targets([NEW_DB]),
        f"{RES}/{PREFIX}_{NEW_DB}_no_identity.tsv"


include: "rules/blast.smk"
include: "rules/annotate.smk"
include: "rules/plot.smk"


# Sites still unmatched after the update. Expected to be empty; anything left
# is explained in <update_db outdir>/<name>.panisa_rescue.tsv.
rule no_identity_report:
    input:
        f"{RES}/{PREFIX}_{NEW_DB}_master.tsv"
    output:
        f"{RES}/{PREFIX}_{NEW_DB}_no_identity.tsv"
    shell:
        """
        awk -F'\t' 'NR==1 || $6 == "No identity"' {input} > {output}
        echo "$(( $(wc -l < {output}) - 1 )) site(s) still 'No identity' against {NEW_DB}" >&2
        """
