
#!/usr/bin/env python3
"""Validate the standalone LDSC results reported in Figure 3.1.

Checks SNP-based heritability, genetic correlations, confidence intervals,
and BH-FDR-adjusted p-values using summary-level results.

Optionally compares estimates against the original LDSC execution logs.
Uses only Python's standard library and does not rerun LDSC.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


# ---------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
OUT = BASE / "results"


# ---------------------------------------------------------------------
# Phenotype definitions
# ---------------------------------------------------------------------

# Seven epilepsy subtypes used in the genetic correlation analyses.
EPI = {
    "CAE",
    "GTCSA",
    "JAE",
    "JME",
    "F-HS",
    "F-Neg",
    "F-Oth",
}

# Seven subcortical brain-volume phenotypes.
BRAIN = {
    "ACC",
    "AMY",
    "BRS",
    "CAU",
    "PAL",
    "PUT",
    "THA",
}

# Phenotype order used in the genetic correlation matrix.
# Broad GGE and FE are excluded from these analyses.
ORDER = [
    "ACC",
    "AMY",
    "BRS",
    "CAU",
    "PAL",
    "PUT",
    "THA",
    "CAE",
    "GTCSA",
    "JAE",
    "JME",
    "F-HS",
    "F-Neg",
    "F-Oth",
]

# Map original LDSC heritability log filenames to figure labels.
# Includes nine epilepsy phenotypes and seven brain volumes.
LOG_NAMES = {
    "GGE": "GGE",
    "FE": "FE",
    "CAE": "CAE",
    "GTCSA": "GTCSA",
    "JAE": "JAE",
    "JME": "JME",
    "FE_HS": "F-HS",
    "FE_Neg": "F-Neg",
    "FE_Oth": "F-Oth",
    "accumbens": "ACC",
    "amygdala": "AMY",
    "brainstem": "BRS",
    "caudate": "CAU",
    "pallidum": "PAL",
    "putamen": "PUT",
    "thalamus": "THA",
}


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def rows(path: Path, delimiter=","):
    """Read a delimited file as a list of dictionaries."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def pair_key(a, b):
    """Return an order-independent phenotype pair."""
    return tuple(sorted((a, b)))


def check(condition, reason):
    """Raise an error when a validation condition is not satisfied."""
    if not condition:
        raise AssertionError(reason)


def approx(a, b, tol=0.00011):
    """Compare numerical values within an absolute tolerance."""
    return math.isclose(
        float(a),
        float(b),
        abs_tol=tol,
        rel_tol=0,
    )


def fdr_bh(pvalues):
    """Calculate Benjamini-Hochberg FDR-adjusted p-values.

    Adjusted p-values are returned in the original input order.
    """
    n = len(pvalues)
    indices = sorted(range(n), key=lambda i: pvalues[i])

    out = [float("nan")] * n
    minimum = 1.0

    for rank_position in range(n - 1, -1, -1):
        idx = indices[rank_position]

        minimum = min(
            minimum,
            pvalues[idx] * n / (rank_position + 1),
            1.0,
        )

        out[idx] = minimum

    return out


# ---------------------------------------------------------------------
# Figure 3.1A: SNP-based heritability
# ---------------------------------------------------------------------

def check_heritability(log_dir=None):
    """Validate heritability estimates and confidence intervals.

    When a log directory is supplied, compare the 16 figure estimates
    against their original standalone LDSC execution logs.
    """

    df = rows(DATA / "figure_3_1A_heritability.csv")

    # Nine epilepsy phenotypes and seven brain volumes.
    check(
        len(df) == 16,
        "Expected 16 heritability entries: 9 epilepsy and 7 volumes",
    )

    name_map = {r["Phenotype"]: r for r in df}

    check(
        len(name_map) == len(df),
        "Duplicate heritability phenotype",
    )

    check(
        set(name_map) == EPI | BRAIN | {"GGE", "FE"},
        "Unexpected heritability phenotypes",
    )

    # Verify 95% confidence intervals.
    for r in df:
        phenotype = r["Phenotype"]
        h2 = float(r["h2"])
        se = float(r["SE"])

        check(
            approx(
                float(r["CI_lower"]),
                h2 - 1.96 * se,
                0.0000001,
            ),
            "Incorrect lower CI for " + phenotype,
        )

        check(
            approx(
                float(r["CI_upper"]),
                h2 + 1.96 * se,
                0.0000001,
            ),
            "Incorrect upper CI for " + phenotype,
        )

    # Confirm broad-phenotype estimates reported in Figure 3.1A.
    check(
        approx(name_map["GGE"]["h2"], 0.2317, 1e-8),
        "Unexpected broad GGE heritability estimate",
    )

    check(
        approx(name_map["GGE"]["SE"], 0.0178, 1e-8),
        "Unexpected broad GGE standard error",
    )

    check(
        approx(name_map["FE"]["h2"], 0.0287, 1e-8),
        "Unexpected broad FE heritability estimate",
    )

    check(
        approx(name_map["FE"]["SE"], 0.0064, 1e-8),
        "Unexpected broad FE standard error",
    )

    # Optionally compare figure values with the original LDSC logs.
    if log_dir is not None:
        count = 0

        for stem, display in LOG_NAMES.items():
            logfile = log_dir / (stem + ".log")

            check(
                logfile.is_file(),
                "Missing heritability log: " + str(logfile),
            )

            content = logfile.read_text(encoding="utf-8")

            match = re.search(
                r"Total (?:Liability|Observed) scale h2:"
                r"\s*([0-9.eE+-]+)\s*\(([0-9.eE+-]+)\)",
                content,
            )

            check(
                match is not None,
                "No heritability estimate found in " + str(logfile),
            )

            check(
                approx(match[1], name_map[display]["h2"]),
                "Heritability disagreement for " + stem,
            )

            check(
                approx(match[2], name_map[display]["SE"]),
                "Standard error disagreement for " + stem,
            )

            count += 1

        check(
            count == 16,
            "Expected 16 heritability logs",
        )

        print(
            f"Validated {count} original standalone LDSC "
            "heritability logs against Figure 3.1A."
        )

    else:
        print(
            "Validated 16 Figure 3.1A heritability estimates "
            "and their confidence intervals. "
            "Original LDSC logs were not supplied."
        )

    return df


# ---------------------------------------------------------------------
# Figure 3.1B: Pairwise genetic correlations
# ---------------------------------------------------------------------

def check_correlations(log_dir=None):
    """Validate pairwise genetic correlations and BH-FDR results."""

    rg = rows(
        DATA / "figure_3_1B_rg_summary.tsv",
        "\t",
    )

    check(
        len(rg) == 91,
        "Expected 91 unique correlations among 14 phenotypes",
    )

    # Validate unique phenotype pairs.
    by_pair = {}

    for r in rg:
        a = r["trait1"]
        b = r["trait2"]

        check(
            a in ORDER and b in ORDER and a != b,
            "Invalid phenotype pair: " + a + ", " + b,
        )

        key = pair_key(a, b)

        check(
            key not in by_pair,
            "Duplicate phenotype pair: " + str(key),
        )

        by_pair[key] = r

    check(
        len(by_pair) == 14 * 13 // 2,
        "Incomplete genetic correlation table",
    )

    # Check consistency between long-format results and the
    # 14 x 14 correlation matrix.
    matrix_path = DATA / "figure_3_1B_rg_matrix.csv"

    with matrix_path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        matrix = list(csv.reader(handle))

    check(
        len(matrix) == 15 and len(matrix[0]) == 15,
        "Expected 14 x 14 genetic correlation matrix",
    )

    check(
        matrix[0][1:] == ORDER,
        "Unexpected column order in genetic correlation matrix",
    )

    for i, row in enumerate(matrix[1:]):
        check(
            row[0] == ORDER[i] and len(row) == 15,
            "Unexpected row in genetic correlation matrix",
        )

        for j, val in enumerate(row[1:]):
            if i == j:
                check(
                    approx(val, 1, 0),
                    "Genetic correlation matrix diagonal is not 1",
                )
            else:
                key = pair_key(ORDER[i], ORDER[j])

                check(
                    approx(
                        val,
                        by_pair[key]["rg"],
                        0.0000001,
                    ),
                    "Matrix versus long-format rg discrepancy",
                )

    # Identify the 49 cross-domain epilepsy-brain correlations.
    cross = [
        r for r in rg
        if (
            r["trait1"] in EPI and r["trait2"] in BRAIN
        ) or (
            r["trait2"] in EPI and r["trait1"] in BRAIN
        )
    ]

    check(
        len(cross) == 49,
        "Expected 49 epilepsy-brain-volume correlations",
    )

    # Recalculate BH-FDR-adjusted p-values.
    qvals = fdr_bh([
        float(r["p"])
        for r in cross
    ])

    # Check nominal associations against the reported results.
    nominal = [
        (r["trait1"], r["trait2"])
        for r in cross
        if float(r["p"]) < 0.05
    ]

    expected_nominal = {
        ("GTCSA", "AMY"),
        ("JAE", "CAU"),
        ("JME", "BRS"),
        ("JME", "PAL"),
    }

    check(
        set(nominal) == expected_nominal,
        "Nominal associations differ from the reported results",
    )

    check(
        all(q >= 0.05 for q in qvals),
        "Unexpected FDR-significant genetic correlation",
    )

    # Verify the adjusted p-values for the four nominal pairs.
    check(
        all(
            approx(q, 0.521850375, 0.00001)
            for r, q in zip(cross, qvals)
            if (r["trait1"], r["trait2"]) in expected_nominal
        ),
        "Unexpected FDR-adjusted p-value for nominal pairs",
    )

    # Export the BH-FDR-adjusted cross-domain results.
    OUT.mkdir(parents=True, exist_ok=True)

    destination = OUT / "figure_3_1_cross_domain_fdr.csv"

    with destination.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.writer(handle)

        writer.writerow([
            "epilepsy_phenotype",
            "brain_volume",
            "rg",
            "se",
            "z",
            "p_unadjusted",
            "q_BH_FDR",
        ])

        for r, q in zip(cross, qvals):
            ep = (
                r["trait1"]
                if r["trait1"] in EPI
                else r["trait2"]
            )

            bv = (
                r["trait2"]
                if r["trait2"] in BRAIN
                else r["trait1"]
            )

            writer.writerow([
                ep,
                bv,
                r["rg"],
                r["se"],
                r["z"],
                r["p"],
                format(q, ".12g"),
            ])

    # Optionally compare all pairwise estimates with original LDSC logs.
    if log_dir is not None:
        observed = {}

        logfiles = sorted(log_dir.glob("*_rg.log"))

        check(
            len(logfiles) == 13,
            "Expected 13 standalone LDSC genetic correlation logs",
        )

        for logfile in logfiles:
            lines = logfile.read_text(
                encoding="utf-8"
            ).splitlines()

            try:
                marker = next(
                    i for i, line in enumerate(lines)
                    if line.startswith(
                        "Summary of Genetic Correlation Results"
                    )
                )
            except StopIteration:
                raise AssertionError(
                    "Missing genetic correlation summary in "
                    + str(logfile)
                )

            for line in lines[marker + 2:]:
                parts = line.split()

                if (
                    len(parts) < 6
                    or not parts[0].endswith(".sumstats.gz")
                    or not parts[1].endswith(".sumstats.gz")
                ):
                    continue

                stems = [
                    Path(path).name.replace(
                        ".hm3.sumstats.gz",
                        "",
                    )
                    for path in parts[:2]
                ]

                # Match munged file stems to phenotype labels.
                labels = []

                for stem in stems:
                    display_name = None

                    for original, display in LOG_NAMES.items():
                        if (
                            stem == original
                            or stem.startswith(original + "_cau_")
                            or stem.startswith(original + "_21_")
                            or (
                                original == "FE_HS"
                                and stem.startswith("focal_HS_")
                            )
                            or (
                                original == "FE_Neg"
                                and stem.startswith(
                                    "focal_lesion_negative_"
                                )
                            )
                            or (
                                original == "FE_Oth"
                                and stem.startswith(
                                    "focal_lesion_other_"
                                )
                            )
                        ):
                            display_name = display
                            break

                    labels.append(display_name)

                a, b = labels

                check(
                    a is not None and b is not None,
                    "Unknown phenotype names: " + str(stems),
                )

                check(
                    a in ORDER and b in ORDER,
                    "Unexpected phenotype in genetic correlation logs",
                )

                key = pair_key(a, b)

                check(
                    key not in observed,
                    "Duplicate correlation in LDSC logs: "
                    + str(key),
                )

                observed[key] = [
                    float(value)
                    for value in parts[2:6]
                ]

        check(
            len(observed) == 91,
            f"Expected 91 LDSC log estimates, "
            f"found {len(observed)}",
        )

        for key, (value, se, z, pval) in observed.items():
            check(
                key in by_pair,
                "Unexpected phenotype pair in LDSC logs: "
                + str(key),
            )

            ref = by_pair[key]

            for label, actual in zip(
                ["rg", "se", "z", "p"],
                [value, se, z, pval],
            ):
                check(
                    approx(
                        actual,
                        ref[label],
                        0.0000001,
                    ),
                    f"{key}: mismatched {label}",
                )

        print(
            "Validated 91 pairwise estimates "
            "(rg/SE/Z/p) against 13 original LDSC logs."
        )

    else:
        print(
            "Validated 91 genetic correlation estimates "
            "against the saved matrix. "
            "Original LDSC logs were not supplied."
        )

    print(
        "BH-FDR: 49 cross-domain tests, "
        "four nominal associations, "
        "zero q < 0.05."
    )

    print(
        "Wrote adjusted results to " + str(destination)
    )


# ---------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--h2-logs",
        type=Path,
        help=(
            "Directory containing the 16 original standalone "
            "LDSC heritability logs"
        ),
    )

    parser.add_argument(
        "--rg-logs",
        type=Path,
        help=(
            "Directory containing the 13 original standalone "
            "LDSC genetic correlation logs"
        ),
    )

    args = parser.parse_args()

    check_heritability(args.h2_logs)
    check_correlations(args.rg_logs)

    print(
        "PASS: Standalone LDSC summary-level results "
        "are internally consistent."
    )


if __name__ == "__main__":
    main()
