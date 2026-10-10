# Step 5: Call IS insertion sites with panISa

rule panISa:
    input:
        final_bam = f"{RES}/{{sample}}/align_reads/{{sample}}_final.bam",
        final_bai = f"{RES}/{{sample}}/align_reads/{{sample}}_final.bam.bai"
    output:
        panISa_out = f"{RES}/{{sample}}/panISa/{{sample}}_panISa.txt"
    params:
        extra = config.get("panisa_args", "")
    log:
        f"{LOGS}/{{sample}}/panISa/{{sample}}_panISa.log"
    conda:
        "../envs/panisa.yaml"
    shell:
        "panISa.py {input.final_bam} -o {output.panISa_out} {params.extra} &> {log}"
