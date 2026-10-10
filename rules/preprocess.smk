# Read preprocessing, same steps as ISScreener:
# trimmomatic -> bwa mem -> picard MarkDuplicates -> samtools sort/index

TRIM = config["trimmomatic"]


def trimmed_reads(wildcards):
    if config.get("trim", True):
        return {"r1": f"{RES}/{wildcards.sample}/trimmomatic/{wildcards.sample}_R1_trim_paired.fastq.gz",
                "r2": f"{RES}/{wildcards.sample}/trimmomatic/{wildcards.sample}_R2_trim_paired.fastq.gz"}
    # pretrimmed mode (ISScreener_pretrimmed.smk)
    return {"r1": config["trimmed_reads"].format(sample=wildcards.sample, read=1),
            "r2": config["trimmed_reads"].format(sample=wildcards.sample, read=2)}


# Step 1: Trim raw fastq files using trimmomatic
rule trim_raw_reads:
    input:
        r1 = os.path.join(config["input_reads"], "{sample}_R1.fastq.gz"),
        r2 = os.path.join(config["input_reads"], "{sample}_R2.fastq.gz")
    output:
        r1_paired = f"{RES}/{{sample}}/trimmomatic/{{sample}}_R1_trim_paired.fastq.gz",
        r2_paired = f"{RES}/{{sample}}/trimmomatic/{{sample}}_R2_trim_paired.fastq.gz",
        r1_unpaired = f"{RES}/{{sample}}/trimmomatic/{{sample}}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = f"{RES}/{{sample}}/trimmomatic/{{sample}}_R2_trim_unpaired.fastq.gz"
    params:
        illuminaclip = "{}:{}:{}:{}:{}:{}".format(
            TRIM["adapter_file_path"], TRIM["seeds"], TRIM["palindrome_clip_threshold"],
            TRIM["simple_clip_threshold"], TRIM["min_adapter_length"],
            str(TRIM["keep_both_reads"]).lower()),
        window = f"{TRIM['window_size']}:{TRIM['window_size_quality']}",
        min_length = TRIM["min_length"],
        head_crop_length = TRIM["head_crop_length"]
    threads: THREADS
    log:
        f"{LOGS}/{{sample}}/trimmomatic/{{sample}}_trimmomatic.log"
    conda:
        "../envs/preprocess.yaml"
    shell:
        "trimmomatic PE -threads {threads} {input.r1} {input.r2} "
        "{output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} "
        "ILLUMINACLIP:{params.illuminaclip} SLIDINGWINDOW:{params.window} "
        "MINLEN:{params.min_length} HEADCROP:{params.head_crop_length} &> {log}"


# Index the reference (skipped when the index is already there)
rule bwa_index:
    input:
        REF
    output:
        multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        f"{LOGS}/reference/bwa_index.log"
    conda:
        "../envs/preprocess.yaml"
    shell:
        "bwa index {input} &> {log}"


rule faidx:
    input:
        REF
    output:
        REF + ".fai"
    conda:
        "../envs/preprocess.yaml"
    shell:
        "samtools faidx {input}"


# Step 2: Align reads using bwa and sort with samtools
rule align_reads:
    input:
        unpack(trimmed_reads),
        ref = REF,
        index = multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        sorted_bam = temp(f"{RES}/{{sample}}/align_reads/{{sample}}_sorted_aln.bam")
    params:
        outdir_temp = f"{RES}/{{sample}}/align_reads/{{sample}}_sorted_aln_temp"
    threads: THREADS
    log:
        f"{LOGS}/{{sample}}/align_reads/{{sample}}_bwa.log"
    conda:
        "../envs/preprocess.yaml"
    shell:
        'bash "{SCRIPTS}/bwa.sh" {input.r1} {input.r2} {input.ref} {threads} '
        "{output.sorted_bam} {params.outdir_temp} {log}"


# Step 3: Remove PCR duplicates
rule remove_pcr_duplicates:
    input:
        sorted_bam = f"{RES}/{{sample}}/align_reads/{{sample}}_sorted_aln.bam"
    output:
        bam_duplicates_removed = temp(f"{RES}/{{sample}}/align_reads/{{sample}}_aln_marked.bam"),
        picard_metrics = f"{RES}/{{sample}}/align_reads/{{sample}}_picard_metrics.txt"
    log:
        f"{LOGS}/{{sample}}/align_reads/{{sample}}_picard.log"
    conda:
        "../envs/preprocess.yaml"
    shell:
        "picard MarkDuplicates -REMOVE_DUPLICATES true -INPUT {input.sorted_bam} "
        "-OUTPUT {output.bam_duplicates_removed} -METRICS_FILE {output.picard_metrics} "
        "-CREATE_INDEX false -VALIDATION_STRINGENCY LENIENT &> {log}"


# Step 4: Sort and index the deduplicated bam
rule bam_sort:
    input:
        bam_duplicates_removed = f"{RES}/{{sample}}/align_reads/{{sample}}_aln_marked.bam"
    output:
        final_bam = f"{RES}/{{sample}}/align_reads/{{sample}}_final.bam",
        final_bai = f"{RES}/{{sample}}/align_reads/{{sample}}_final.bam.bai"
    params:
        samtools_temp = f"{RES}/{{sample}}/align_reads/{{sample}}_final_temp"
    threads: THREADS
    log:
        f"{LOGS}/{{sample}}/align_reads/{{sample}}_post_picard_sort.log"
    conda:
        "../envs/preprocess.yaml"
    shell:
        """
        samtools sort {input.bam_duplicates_removed} -m 500M -@ {threads} -o {output.final_bam} -T {params.samtools_temp} &> {log}
        samtools index {output.final_bam}
        samtools flagstat {output.final_bam} >> {log}
        """
