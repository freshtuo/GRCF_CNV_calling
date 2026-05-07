# Optional FACETS tumor-normal analysis for human paired comparisons.
# Eligibility is enforced by metadata validation before this rule can run.
rule facets:
    input:
        validated=f"{RESULTS}/metadata/validated.ok"
    output:
        segments=f"{RESULTS}/facets/{{comparison_id}}/segments/facets_segments.tsv",
        purity_ploidy=f"{RESULTS}/facets/{{comparison_id}}/purity_ploidy/facets_purity_ploidy.tsv",
        plot=f"{RESULTS}/facets/{{comparison_id}}/plots/facets.pdf"
    params:
        common_snps=lambda wildcards: config["resources"]["human"].get("common_snps_vcf", ""),
        case_bam=lambda wildcards: sample_bam(comparison_case_id(wildcards.comparison_id)),
        control_bam=lambda wildcards: sample_bam(comparison_control_id(wildcards.comparison_id)),
        outdir=lambda wildcards, output: str(Path(output.segments).parents[1]),
        comparison_id="{comparison_id}"
    log:
        f"{RESULTS}/facets/{{comparison_id}}/facets.log"
    conda:
        "envs/facets.yaml"
    shell:
        r"""
        set -euo pipefail
        # FACETS requires a common SNP VCF to build tumor/normal allele counts.
        test -n "{params.common_snps}"
        mkdir -p {params.outdir}/segments {params.outdir}/purity_ploidy {params.outdir}/plots {params.outdir}/pileup
        # snp-pileup creates the allele-count matrix consumed by FACETS in R.
        snp-pileup --gzip --pseudo-snps 100 {params.common_snps} \
            {params.outdir}/pileup/{params.comparison_id}.csv.gz \
            {params.control_bam} {params.case_bam} > {log} 2>&1
        # The R helper writes segment calls, purity/ploidy, and a diagnostic plot.
        Rscript scripts/run_facets.R \
            --pileup {params.outdir}/pileup/{params.comparison_id}.csv.gz \
            --comparison-id {params.comparison_id} \
            --segments {output.segments} \
            --purity-ploidy {output.purity_ploidy} \
            --plot {output.plot} >> {log} 2>&1
        """
