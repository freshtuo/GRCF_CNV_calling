import csv
from pathlib import Path


# Shared values and helper functions used by all modular rule files.
# This file is included first by Snakefile so later rules can ask questions like
# "which BAM belongs to this comparison?" without duplicating metadata parsing.
PROJECT = config["project"]
OUTDIR = config.get("outdir", "results")
RESULTS = f"{OUTDIR}/{PROJECT}"
SAMPLES_TSV = config.get("samples_tsv", "config/samples.tsv")
COMPARISONS_TSV = config.get("comparisons_tsv", "config/comparisons.tsv")


def _read_tsv(path):
    """Load a tab-delimited metadata table as a list of dictionaries."""
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _is_yes(value):
    """Interpret yes/true/1/y values from config tables as enabled flags."""
    return str(value).strip().lower() in {"yes", "true", "1", "y"}


# samples.tsv describes physical samples; comparisons.tsv decides analyses.
SAMPLES = {row["sample_id"]: row for row in _read_tsv(SAMPLES_TSV)}
COMPARISONS = {row["comparison_id"]: row for row in _read_tsv(COMPARISONS_TSV)}

# These lists determine which comparison IDs appear in CNVkit and FACETS targets.
CNVKIT_COMPARISONS = [
    cid for cid, row in COMPARISONS.items() if _is_yes(row.get("run_cnvkit", ""))
]
FACETS_COMPARISONS = [
    cid for cid, row in COMPARISONS.items() if _is_yes(row.get("run_facets", ""))
]


def sample_bam(sample_id):
    """Return the BAM path for a sample_id from samples.tsv."""
    return SAMPLES[sample_id]["bam"]


def sample_bai(sample_id):
    """Return the BAI path for a sample_id from samples.tsv."""
    return SAMPLES[sample_id]["bai"]


def comparison_case_id(comparison_id):
    """Return the case sample_id for a comparison_id."""
    return COMPARISONS[comparison_id]["case_id"]


def comparison_control_id(comparison_id):
    """Return the optional control sample_id for a comparison_id."""
    return COMPARISONS[comparison_id].get("control_id", "").strip()


def comparison_species(comparison_id):
    """Return the species label for a comparison_id."""
    return COMPARISONS[comparison_id]["species"]


def species_resource(comparison_id, key):
    """Return the configured reference resource for a comparison's species."""
    species = comparison_species(comparison_id)
    return config["resources"][species][key]


# Standard CNVkit outputs are named by comparison_id so each comparison folder is
# self-contained even when the same case sample is analyzed against several controls.
CNVKIT_CNR = f"{RESULTS}/cnvkit/{{comparison_id}}/cnr/{{comparison_id}}.cnr"
CNVKIT_CNS = f"{RESULTS}/cnvkit/{{comparison_id}}/cns/{{comparison_id}}.cns"
CNVKIT_CALL_CNS = f"{RESULTS}/cnvkit/{{comparison_id}}/calls/{{comparison_id}}.call.cns"
