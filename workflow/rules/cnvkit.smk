# Run CNVkit coverage/binning/segmentation for one comparison.
# If control_id exists, normal_arg adds `--normal`; otherwise CNVkit runs tumor-only.
rule cnvkit_batch:
    input:
        validated=f"{RESULTS}/metadata/validated.ok"
    output:
        cnr=CNVKIT_CNR,
        cns=CNVKIT_CNS
    threads:
        int(config.get("cnvkit", {}).get("threads", 1))
    resources:
        mem_mb=int(config.get("cnvkit", {}).get("mem_mb", 8000))
    params:
        case_id=lambda wildcards: comparison_case_id(wildcards.comparison_id),
        case_bam=lambda wildcards: sample_bam(comparison_case_id(wildcards.comparison_id)),
        normal_arg=lambda wildcards: f"--normal {sample_bam(comparison_control_id(wildcards.comparison_id))}" if comparison_control_id(wildcards.comparison_id) else "",
        fasta=lambda wildcards: species_resource(wildcards.comparison_id, "fasta"),
        access=lambda wildcards: species_resource(wildcards.comparison_id, "access_bed"),
        annotate=lambda wildcards: species_resource(wildcards.comparison_id, "cnvkit_annotate"),
        outdir=lambda wildcards, output: str(Path(output.cnr).parents[1]),
        batch_dir=lambda wildcards, output: str(Path(output.cnr).parents[1] / "batch_outputs"),
        case_prefix=lambda wildcards: Path(sample_bam(comparison_case_id(wildcards.comparison_id))).stem
    log:
        f"{RESULTS}/cnvkit/{{comparison_id}}/cnvkit_batch.log"
    conda:
        "envs/cnvkit.yaml"
    shell:
        r"""
        set -euo pipefail
        # Keep the public result layout stable while allowing CNVkit to write
        # its native filenames inside a temporary working directory.
        rm -rf {params.batch_dir}
        mkdir -p {params.outdir}/cnr {params.outdir}/cns {params.outdir}/calls {params.outdir}/plots {params.outdir}/gene_annotation {params.batch_dir}
        workdir=$(mktemp -d {params.outdir}/batch_tmp.XXXXXX)
        cnvkit.py batch {params.case_bam} {params.normal_arg} \
            --method wgs \
            --fasta {params.fasta} \
            --access {params.access} \
            --annotate {params.annotate} \
            --processes {threads} \
            --output-dir "$workdir" > {log} 2>&1
        # CNVkit names outputs from the input BAM. Preserve the native output
        # set, then copy canonical comparison_id names for downstream rules.
        cnr="$workdir/{params.case_prefix}.cnr"
        cns="$workdir/{params.case_prefix}.cns"
        test -f "$cnr"
        test -f "$cns"
        cp -a "$workdir"/. {params.batch_dir}/
        cp "$cnr" {output.cnr}
        cp "$cns" {output.cns}
        rm -rf "$workdir"
        """


# Convert segmented log2 ratios into called copy-number states.
rule cnvkit_call:
    input:
        cns=CNVKIT_CNS
    output:
        call_cns=CNVKIT_CALL_CNS
    resources:
        mem_mb=int(config.get("cnvkit", {}).get("mem_mb", 8000))
    log:
        f"{RESULTS}/cnvkit/{{comparison_id}}/calls/cnvkit_call.log"
    conda:
        "envs/cnvkit.yaml"
    shell:
        "cnvkit.py call {input.cns} -o {output.call_cns} > {log} 2>&1"


# Produce a genome-wide scatter plot for visual review of bins and segments.
rule cnvkit_scatter:
    input:
        cnr=CNVKIT_CNR,
        cns=CNVKIT_CNS
    output:
        f"{RESULTS}/cnvkit/{{comparison_id}}/plots/{'{comparison_id}'}.scatter.pdf"
    resources:
        mem_mb=int(config.get("cnvkit", {}).get("mem_mb", 8000))
    log:
        f"{RESULTS}/cnvkit/{{comparison_id}}/plots/scatter.log"
    conda:
        "envs/cnvkit.yaml"
    shell:
        "cnvkit.py scatter {input.cnr} -s {input.cns} -o {output} > {log} 2>&1"


# Produce a chromosome-level CNV diagram from the CNVkit segments.
rule cnvkit_diagram:
    input:
        cns=CNVKIT_CNS
    output:
        f"{RESULTS}/cnvkit/{{comparison_id}}/plots/{'{comparison_id}'}.diagram.pdf"
    resources:
        mem_mb=int(config.get("cnvkit", {}).get("mem_mb", 8000))
    log:
        f"{RESULTS}/cnvkit/{{comparison_id}}/plots/diagram.log"
    conda:
        "envs/cnvkit.yaml"
    shell:
        "cnvkit.py diagram {input.cns} -o {output} > {log} 2>&1"
