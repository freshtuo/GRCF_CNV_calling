# Annotate CNVkit called segments by overlapping them with the species-specific
# genes_bed and applying CNV/event-size thresholds from config.yaml.
rule annotate_cnvkit:
    input:
        call_cns=CNVKIT_CALL_CNS,
        cns=CNVKIT_CNS,
        scatter=f"{RESULTS}/cnvkit/{{comparison_id}}/plots/{'{comparison_id}'}.scatter.png",
        diagram=f"{RESULTS}/cnvkit/{{comparison_id}}/plots/{'{comparison_id}'}.diagram.pdf"
    output:
        segments=CNVKIT_ANNOTATED_SEGMENTS,
        genes=CNVKIT_ANNOTATED_GENES
    params:
        genes_bed=lambda wildcards: species_resource(wildcards.comparison_id, "genes_bed"),
        thresholds=lambda wildcards: ",".join([
            str(config["cnv_thresholds"]["gain"]),
            str(config["cnv_thresholds"]["loss"]),
            str(config["cnv_thresholds"]["amp"]),
            str(config["cnv_thresholds"]["deep_del"]),
            str(config["event_size"]["focal_mb"]),
            str(config["event_size"]["broad_mb"]),
        ]),
        min_gene_overlap_bp=lambda wildcards: config.get("annotation", {}).get("min_gene_overlap_bp", 1),
        min_gene_overlap_fraction=lambda wildcards: config.get("annotation", {}).get("min_gene_overlap_fraction", 0.0)
    log:
        f"{RESULTS}/cnvkit/{{comparison_id}}/annotation/annotation.log"
    conda:
        "../../envs/annotation.yaml"
    shell:
        "mkdir -p $(dirname {log}); "
        "python scripts/annotate_segments.py "
        "--source-tool cnvkit "
        "--segments {input.call_cns} "
        "--genes-bed {params.genes_bed} "
        "--thresholds {params.thresholds} "
        "--min-gene-overlap-bp {params.min_gene_overlap_bp} "
        "--min-gene-overlap-fraction {params.min_gene_overlap_fraction} "
        "--out-segments {output.segments} "
        "--out-genes {output.genes} > {log} 2>&1"


# Annotate FACETS segments when optional FACETS analysis is enabled.
# The same annotation script normalizes FACETS columns into the shared schema.
rule annotate_facets:
    input:
        segments=FACETS_SEGMENTS
    output:
        segments=FACETS_ANNOTATED_SEGMENTS,
        genes=FACETS_ANNOTATED_GENES
    params:
        genes_bed=lambda wildcards: species_resource(wildcards.comparison_id, "genes_bed"),
        thresholds=lambda wildcards: ",".join([
            str(config["cnv_thresholds"]["gain"]),
            str(config["cnv_thresholds"]["loss"]),
            str(config["cnv_thresholds"]["amp"]),
            str(config["cnv_thresholds"]["deep_del"]),
            str(config["event_size"]["focal_mb"]),
            str(config["event_size"]["broad_mb"]),
        ]),
        min_gene_overlap_bp=lambda wildcards: config.get("annotation", {}).get("min_gene_overlap_bp", 1),
        min_gene_overlap_fraction=lambda wildcards: config.get("annotation", {}).get("min_gene_overlap_fraction", 0.0)
    log:
        f"{RESULTS}/facets/{{comparison_id}}/annotation/annotation.log"
    conda:
        "../../envs/annotation.yaml"
    shell:
        "mkdir -p $(dirname {log}); "
        "python scripts/annotate_segments.py "
        "--source-tool facets "
        "--segments {input.segments} "
        "--genes-bed {params.genes_bed} "
        "--thresholds {params.thresholds} "
        "--min-gene-overlap-bp {params.min_gene_overlap_bp} "
        "--min-gene-overlap-fraction {params.min_gene_overlap_fraction} "
        "--out-segments {output.segments} "
        "--out-genes {output.genes} > {log} 2>&1"
