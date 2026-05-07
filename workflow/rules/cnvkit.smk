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
        "../../envs/cnvkit.yaml"
    shell:
        r"""
        set -euo pipefail
        # Keep the public result layout stable while allowing CNVkit to write
        # its native filenames inside a temporary working directory.
        rm -rf {params.batch_dir}
        mkdir -p {params.outdir}/cnr {params.outdir}/cns {params.outdir}/calls {params.outdir}/plots {params.outdir}/annotation {params.batch_dir}
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
        chmod -R u+rwX,g+rX {params.outdir}
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
        "../../envs/cnvkit.yaml"
    shell:
        "cnvkit.py call {input.cns} -o {output.call_cns} > {log} 2>&1"


# Produce a genome-wide scatter plot for visual review of bins and segments.
rule cnvkit_scatter:
    input:
        cnr=CNVKIT_CNR,
        cns=CNVKIT_CNS
    output:
        pdf=f"{RESULTS}/cnvkit/{{comparison_id}}/plots/{{comparison_id}}.scatter.pdf",
        png=f"{RESULTS}/cnvkit/{{comparison_id}}/plots/{{comparison_id}}.scatter.png",
        chrom_done=f"{RESULTS}/cnvkit/{{comparison_id}}/plots/chromosomes/{{comparison_id}}.done"
    resources:
        mem_mb=int(config.get("cnvkit", {}).get("mem_mb", 8000))
    log:
        f"{RESULTS}/cnvkit/{{comparison_id}}/plots/scatter.log"
    conda:
        "../../envs/cnvkit.yaml"
    params:
        fig_size=config.get("cnvkit", {}).get("scatter_fig_size", "18 5"),
        chrom_fig_size=config.get("cnvkit", {}).get("chromosome_scatter_fig_size", "14 5")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.pdf})" "$(dirname {output.chrom_done})"
        cnvkit.py scatter {input.cnr} -s {input.cns} --fig-size {params.fig_size} -o {output.pdf} > {log} 2>&1
        cnvkit.py scatter {input.cnr} -s {input.cns} --fig-size {params.fig_size} -o {output.png} >> {log} 2>&1
        cut -f1 {input.cnr} | tail -n +2 | sort -u | while read -r chrom; do
            test -n "$chrom" || continue
            cnvkit.py scatter {input.cnr} -s {input.cns} -c "$chrom" --fig-size {params.chrom_fig_size} \
                -o "$(dirname {output.chrom_done})/{wildcards.comparison_id}.${{chrom}}.scatter.png" >> {log} 2>&1
        done
        touch {output.chrom_done}
        """


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
        "../../envs/cnvkit.yaml"
    shell:
        "cnvkit.py diagram {input.cns} -o {output} > {log} 2>&1"
