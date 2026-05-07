# Validate all user-provided metadata before any BAM, CNVkit, or FACETS work begins.
# The touched .ok file is a small checkpoint that downstream rules depend on.
rule validate_metadata:
    input:
        samples=SAMPLES_TSV,
        comparisons=COMPARISONS_TSV
    output:
        touch(f"{RESULTS}/metadata/validated.ok")
    log:
        f"{RESULTS}/metadata/validation.log"
    params:
        config_yaml="config/config.yaml"
    conda:
        "../../envs/annotation.yaml"
    shell:
        "python scripts/validate_metadata.py "
        "--config {params.config_yaml} "
        "--samples {input.samples} "
        "--comparisons {input.comparisons} "
        "--log {log}"
