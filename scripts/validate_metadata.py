#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
import yaml


REQUIRED_SAMPLE_COLUMNS = {"sample_id", "patient_id", "species", "sample_type", "bam", "bai"}
REQUIRED_COMPARISON_COLUMNS = {
    "comparison_id",
    "patient_id",
    "species",
    "case_id",
    "control_id",
    "comparison_type",
    "run_cnvkit",
    "run_facets",
}


def read_tsv(path):
    """Read a tab-delimited metadata file as strings and keep blanks as empty."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def yes(value):
    """Normalize yes/no style metadata fields to a boolean decision."""
    return str(value).strip().lower() in {"yes", "true", "1", "y"}


def require_columns(kind, actual, required, errors):
    """Record a validation error if a metadata file is missing required columns."""
    missing = sorted(required - actual)
    if missing:
        errors.append(f"{kind} is missing required columns: {', '.join(missing)}")


def main():
    """Validate config, samples.tsv, and comparisons.tsv before workflow jobs run."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--comparisons", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    with open(args.config) as handle:
        cfg = yaml.safe_load(handle)

    samples = read_tsv(args.samples)
    comparisons = read_tsv(args.comparisons)

    errors = []
    warnings = []

    # Check file structure first. Later checks assume these columns are present.
    require_columns("samples.tsv", set(samples.columns), REQUIRED_SAMPLE_COLUMNS, errors)
    require_columns("comparisons.tsv", set(comparisons.columns), REQUIRED_COMPARISON_COLUMNS, errors)
    if errors:
        write_log(args.log, warnings, errors)
        raise SystemExit(f"Metadata validation failed with {len(errors)} error(s); see {args.log}")

    # Build a lookup table so each comparison can resolve case/control sample IDs quickly.
    sample_by_id = samples.set_index("sample_id", drop=False).to_dict(orient="index")

    # Every sample row must point to existing BAM and BAI files because the pipeline
    # starts from already-aligned data and never realigns reads.
    for sample in samples.to_dict(orient="records"):
        sample_id = sample["sample_id"]
        for key in ("bam", "bai"):
            value = sample[key]
            if not value:
                errors.append(f"Sample {sample_id} has an empty {key} path")
            elif not Path(value).exists():
                errors.append(f"Sample {sample_id} {key} does not exist: {value}")

    # comparisons.tsv is the workflow driver. Each row defines one CNV analysis,
    # so validation is done at the comparison level rather than by sample_type.
    for row in comparisons.to_dict(orient="records"):
        cid = row["comparison_id"]
        species = row["species"]
        case_id = row["case_id"]
        control_id = row["control_id"].strip()
        comparison_type = row["comparison_type"]

        # The case sample is mandatory. If it is missing, the rest of this
        # comparison cannot be interpreted safely.
        if case_id not in sample_by_id:
            errors.append(f"Comparison {cid} case_id is not in samples.tsv: {case_id}")
            continue

        # Species must match across the comparison row and all participating samples
        # so the correct reference FASTA/access/gene resources are used.
        if sample_by_id[case_id].get("species") != species:
            errors.append(
                f"Comparison {cid} species is {species}, but case sample {case_id} "
                f"is {sample_by_id[case_id].get('species')}"
            )

        # control_id is optional for tumor-only CNVkit, but if provided it must
        # exist and have the same species as the comparison.
        if control_id:
            if control_id not in sample_by_id:
                errors.append(f"Comparison {cid} control_id is not in samples.tsv: {control_id}")
            elif sample_by_id[control_id].get("species") != species:
                errors.append(
                    f"Comparison {cid} species is {species}, but control sample {control_id} "
                    f"is {sample_by_id[control_id].get('species')}"
                )

        # FACETS is intentionally limited to human paired comparisons in this
        # initial implementation. Mouse and tumor-only analyses stay CNVkit-only.
        if yes(row["run_facets"]):
            if species != "human":
                errors.append(f"Comparison {cid} has run_facets=yes but species is not human")
            if not control_id:
                errors.append(f"Comparison {cid} has run_facets=yes but control_id is empty")
            if comparison_type == "tumor_only":
                errors.append(f"Comparison {cid} has run_facets=yes but comparison_type is tumor_only")
            if not cfg.get("tools", {}).get("facets", False):
                errors.append(f"Comparison {cid} has run_facets=yes but tools.facets is false")
            if not cfg.get("resources", {}).get("human", {}).get("common_snps_vcf"):
                warnings.append(
                    f"Comparison {cid} requires resources.human.common_snps_vcf before FACETS can run"
                )
            elif not Path(cfg["resources"]["human"]["common_snps_vcf"]).exists():
                errors.append(
                    "Comparison "
                    f"{cid} resources.human.common_snps_vcf does not exist: "
                    f"{cfg['resources']['human']['common_snps_vcf']}"
                )

        # The comparison-level run flags are honored only if the corresponding
        # tool is enabled globally in config.yaml.
        if yes(row["run_cnvkit"]) and not cfg.get("tools", {}).get("cnvkit", False):
            errors.append(f"Comparison {cid} has run_cnvkit=yes but tools.cnvkit is false")

        # Each species used by a comparison needs reference resources for CNVkit
        # and gene-overlap annotation.
        if species not in cfg.get("resources", {}):
            errors.append(f"Comparison {cid} species has no resources entry in config.yaml: {species}")
        else:
            species_resources = cfg.get("resources", {}).get(species, {})
            resource_keys = ["fasta", "access_bed", "genes_bed"]
            if yes(row["run_cnvkit"]):
                resource_keys.append("cnvkit_annotate")
            for key in resource_keys:
                value = species_resources.get(key)
                if not value:
                    errors.append(f"Comparison {cid} species {species} is missing resources.{species}.{key}")
                elif not Path(value).exists():
                    errors.append(
                        f"Comparison {cid} resources.{species}.{key} does not exist: {value}"
                    )

    write_log(args.log, warnings, errors)

    # Non-zero exit tells Snakemake that downstream jobs must not run.
    if errors:
        raise SystemExit(f"Metadata validation failed with {len(errors)} error(s); see {args.log}")


def write_log(path, warnings, errors):
    """Write validation warnings/errors to the Snakemake log file."""
    # Write one log file that Snakemake can keep with the project outputs. Warnings
    # are informative; errors stop the workflow before expensive jobs start.
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as log:
        if warnings:
            log.write("WARNINGS\n")
            for warning in warnings:
                log.write(f"- {warning}\n")
            log.write("\n")
        if errors:
            log.write("ERRORS\n")
            for error in errors:
                log.write(f"- {error}\n")
        else:
            log.write("Metadata validation passed.\n")


if __name__ == "__main__":
    main()
