# Optional post-analysis maintenance targets. These are intentionally excluded
# from rule all so compression and cleanup only run when requested explicitly.


rule gzip_cnvkit_batch_outputs:
    input:
        expand(f"{RESULTS}/maintenance/cnvkit/{{comparison_id}}/batch_outputs.gzip.done", comparison_id=CNVKIT_COMPARISONS)


rule gzip_cnvkit_batch_output:
    input:
        report=f"{RESULTS}/summary/report.html",
        cnr=CNVKIT_CNR,
        cns=CNVKIT_CNS
    output:
        touch(f"{RESULTS}/maintenance/cnvkit/{{comparison_id}}/batch_outputs.gzip.done")
    params:
        batch_dir=lambda wildcards, input: str(Path(input.cnr).parents[1] / "batch_outputs")
    log:
        f"{RESULTS}/maintenance/cnvkit/{{comparison_id}}/batch_outputs.gzip.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        exec > "{log}" 2>&1
        if [ -d "{params.batch_dir}" ]; then
            find "{params.batch_dir}" -type f ! -name '*.gz' \
                \( -name '*.cnr' -o -name '*.cns' -o -name '*.cnn' -o -name '*.bed' -o -name '*.tsv' -o -name '*.txt' \) \
                -print0 | xargs -0 -r gzip -f
        else
            echo "No CNVkit batch output directory found: {params.batch_dir}"
        fi
        """


rule clean_cnvkit_batch_outputs:
    input:
        expand(f"{RESULTS}/maintenance/cnvkit/{{comparison_id}}/batch_outputs.clean.done", comparison_id=CNVKIT_COMPARISONS)


rule clean_cnvkit_batch_output:
    input:
        report=f"{RESULTS}/summary/report.html",
        cnr=CNVKIT_CNR,
        cns=CNVKIT_CNS
    output:
        touch(f"{RESULTS}/maintenance/cnvkit/{{comparison_id}}/batch_outputs.clean.done")
    params:
        batch_dir=lambda wildcards, input: str(Path(input.cnr).parents[1] / "batch_outputs")
    log:
        f"{RESULTS}/maintenance/cnvkit/{{comparison_id}}/batch_outputs.clean.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        exec > "{log}" 2>&1
        if [ -d "{params.batch_dir}" ]; then
            rm -rf "{params.batch_dir}"
            echo "Removed CNVkit batch output directory: {params.batch_dir}"
        else
            echo "No CNVkit batch output directory found: {params.batch_dir}"
        fi
        """
