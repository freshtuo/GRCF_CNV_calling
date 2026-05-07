# Confirm each BAM is readable and structurally valid.
# quickcheck writes nothing on success, so the output file is empty when the BAM passes.
rule samtools_quickcheck:
    input:
        validated=f"{RESULTS}/metadata/validated.ok"
    params:
        bam=lambda wildcards: sample_bam(wildcards.sample_id)
    output:
        f"{RESULTS}/qc/samples/{{sample_id}}/quickcheck.txt"
    log:
        f"{RESULTS}/qc/samples/{{sample_id}}/quickcheck.log"
    conda:
        "envs/annotation.yaml"
    shell:
        "mkdir -p $(dirname {output}) $(dirname {log}); "
        "samtools quickcheck -v {params.bam} > {output} 2> {log}"


# Summarize alignment-level counts such as total, mapped, and properly paired reads.
rule samtools_flagstat:
    input:
        validated=f"{RESULTS}/metadata/validated.ok"
    params:
        bam=lambda wildcards: sample_bam(wildcards.sample_id)
    output:
        f"{RESULTS}/qc/samples/{{sample_id}}/flagstat.txt"
    log:
        f"{RESULTS}/qc/samples/{{sample_id}}/flagstat.log"
    conda:
        "envs/annotation.yaml"
    shell:
        "mkdir -p $(dirname {output}) $(dirname {log}); "
        "samtools flagstat {params.bam} > {output} 2> {log}"


# Count mapped/unmapped reads by reference sequence to catch index or contig issues.
rule samtools_idxstats:
    input:
        validated=f"{RESULTS}/metadata/validated.ok"
    params:
        bam=lambda wildcards: sample_bam(wildcards.sample_id)
    output:
        f"{RESULTS}/qc/samples/{{sample_id}}/idxstats.txt"
    log:
        f"{RESULTS}/qc/samples/{{sample_id}}/idxstats.log"
    conda:
        "envs/annotation.yaml"
    shell:
        "mkdir -p $(dirname {output}) $(dirname {log}); "
        "samtools idxstats {params.bam} > {output} 2> {log}"
