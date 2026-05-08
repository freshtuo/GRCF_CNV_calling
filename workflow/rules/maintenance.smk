# Optional post-analysis maintenance targets. These are intentionally excluded
# from rule all so compression and cleanup only run when requested explicitly.


rule standard_cleanup:
    input:
        f"{RESULTS}/maintenance/facets/pileups.clean.done",
        f"{RESULTS}/maintenance/cnvkit/batch_tmp.clean.done",
        expand(f"{RESULTS}/maintenance/cnvkit/{{comparison_id}}/batch_outputs.gzip.done", comparison_id=CNVKIT_COMPARISONS)


rule cleanup_intermediates:
    input:
        f"{RESULTS}/maintenance/facets/pileups.clean.done",
        f"{RESULTS}/maintenance/cnvkit/batch_tmp.clean.done"


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
    threads:
        config.get("maintenance", {}).get("gzip_threads", 4)
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
            if command -v pigz >/dev/null 2>&1; then
                compressor="pigz -p {threads} -f"
            else
                compressor="gzip -f"
                echo "pigz not found; falling back to single-threaded gzip"
            fi
            find "{params.batch_dir}" -type f ! -name '*.gz' \
                \( -name '*.cnr' -o -name '*.cns' -o -name '*.cnn' -o -name '*.bed' -o -name '*.tsv' -o -name '*.txt' \) \
                -print0 | while IFS= read -r -d '' file; do
                    $compressor "$file"
                done
        else
            echo "No CNVkit batch output directory found: {params.batch_dir}"
        fi
        """


rule clean_facets_pileups:
    input:
        report=f"{RESULTS}/summary/report.html"
    output:
        touch(f"{RESULTS}/maintenance/facets/pileups.clean.done")
    params:
        facets_dir=f"{RESULTS}/facets"
    log:
        f"{RESULTS}/maintenance/facets/pileups.clean.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        exec > "{log}" 2>&1
        if [ -d "{params.facets_dir}" ]; then
            find "{params.facets_dir}" -mindepth 2 -maxdepth 2 -type d -name pileup -print -exec rm -rf {{}} +
        else
            echo "No FACETS directory found: {params.facets_dir}"
        fi
        """


rule clean_cnvkit_batch_tmp:
    input:
        report=f"{RESULTS}/summary/report.html"
    output:
        touch(f"{RESULTS}/maintenance/cnvkit/batch_tmp.clean.done")
    params:
        cnvkit_dir=f"{RESULTS}/cnvkit"
    log:
        f"{RESULTS}/maintenance/cnvkit/batch_tmp.clean.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        exec > "{log}" 2>&1
        if [ -d "{params.cnvkit_dir}" ]; then
            find "{params.cnvkit_dir}" -mindepth 2 -maxdepth 2 -type d -name 'batch_tmp.*' -print -exec rm -rf {{}} +
        else
            echo "No CNVkit directory found: {params.cnvkit_dir}"
        fi
        """


rule deep_cleanup:
    input:
        f"{RESULTS}/maintenance/facets/pileups.clean.done",
        f"{RESULTS}/maintenance/cnvkit/batch_tmp.clean.done",
        expand(f"{RESULTS}/maintenance/cnvkit/{{comparison_id}}/batch_outputs.clean.done", comparison_id=CNVKIT_COMPARISONS)


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
