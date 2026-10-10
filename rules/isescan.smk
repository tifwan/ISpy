# Extract IS elements from each sample assembly with ISEScan.
# Output feeds update_db.smk.

rule isescan:
    input:
        assembly = lambda wildcards: config["assembly"].format(sample=wildcards.sample)
    output:
        is_fna = f"{RES}/isescan/{{sample}}/{{sample}}.is.fna"
    params:
        outdir = f"{RES}/isescan/{{sample}}/isescan_out"
    threads: THREADS
    log:
        f"{LOGS}/{{sample}}/isescan/{{sample}}_isescan.log"
    conda:
        "../envs/isescan.yaml"
    shell:
        """
        rm -rf {params.outdir}
        isescan.py --seqfile {input.assembly} --output {params.outdir} --nthread {threads} &> {log}
        # ISEScan nests its output under the input path; collect the IS sequences.
        # No file means no IS elements were found.
        is_fna=$(find {params.outdir} -name '*.is.fna' | head -n 1)
        if [ -n "$is_fna" ]; then cp "$is_fna" {output.is_fna}; else : > {output.is_fna}; fi
        """
