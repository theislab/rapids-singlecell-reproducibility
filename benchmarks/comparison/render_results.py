from __future__ import annotations

import argparse
import csv
import json
from datetime import UTC, datetime
from pathlib import Path

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render reviewer-facing CPU/GPU equivalence results.")
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--execution", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def status(value: bool) -> str:
    return "PASS" if value else "FAIL"


def is_gating(metric: dict) -> bool:
    return bool(metric.get("gating", True))


def gating_metrics(record: dict) -> list[dict]:
    return [metric for metric in record["metrics"] if is_gating(metric)]


def render_csv(summary: dict, output: Path) -> None:
    fields = [
        "method",
        "reference_package",
        "dataset",
        "tier",
        "metric",
        "observed",
        "comparison",
        "tolerance",
        "gating",
        "basis",
        "passed",
        "versions",
    ]
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for record in summary["records"]:
            for metric in record["metrics"]:
                writer.writerow(
                    {
                        "method": record["method"],
                        "reference_package": record["reference_package"],
                        "dataset": record["dataset"],
                        "tier": record["tier"],
                        "metric": metric["metric"],
                        "observed": metric["observed"],
                        "comparison": metric["comparison"],
                        "tolerance": metric["tolerance"],
                        "gating": is_gating(metric),
                        "basis": metric.get("basis", ""),
                        "passed": metric["passed"],
                        "versions": json.dumps(record.get("versions", {}), sort_keys=True),
                    }
                )


def render_plot(summary: dict, output: Path) -> None:
    records = summary["records"]
    labels = [record["method"] for record in records]
    rates = [
        sum(metric["passed"] for metric in gating_metrics(record)) / max(1, len(gating_metrics(record)))
        for record in records
    ]
    colors = ["#238636" if record["passed"] else "#da3633" for record in records]
    height = max(5, 0.38 * len(records))
    figure, axis = plt.subplots(figsize=(10, height), constrained_layout=True)
    axis.barh(labels, rates, color=colors)
    axis.set_xlim(0, 1.02)
    axis.set_xlabel("Fraction of equivalence metrics passing")
    axis.set_title("CPU/GPU equivalence by comparison group")
    axis.axvline(1, color="#57606a", linewidth=0.8)
    axis.invert_yaxis()
    for row, rate in enumerate(rates):
        axis.text(min(rate + 0.01, 0.96), row, f"{rate:.0%}", va="center", fontsize=8)
    figure.savefig(output, dpi=200)
    plt.close(figure)


def render_markdown(summary: dict, execution: dict, output: Path) -> None:
    methods = summary["n_methods"]
    passed_methods = summary.get("n_passed_methods", sum(record["passed"] for record in summary["records"]))
    metrics = summary["n_metrics"]
    passed_metrics = summary.get(
        "n_passed_metrics",
        sum(metric["passed"] for record in summary["records"] for metric in gating_metrics(record)),
    )
    informational = summary.get(
        "n_informational_metrics",
        sum(1 for record in summary["records"] for metric in record["metrics"] if not is_gating(metric)),
    )
    versions = {}
    for record in summary["records"]:
        for package, package_version in record.get("versions", {}).items():
            versions.setdefault(package, set()).add(str(package_version))
    lines = [
        "# CPU/GPU equivalence report",
        "",
        f"Generated {datetime.now(UTC).isoformat()} from isolated comparison processes.",
        "",
        "## Outcome",
        "",
        *(
            [
                f"> Partial run: only {execution.get('n_selected')} of {execution.get('n_available')} comparisons"
                " were executed. Records for the others, if present, come from an earlier run.",
                "",
            ]
            if execution.get("partial")
            else []
        ),
        f"- Overall: **{status(summary['passed'])}**",
        f"- Method groups passing: **{passed_methods}/{methods}**",
        f"- Gating metrics passing: **{passed_metrics}/{metrics}**",
        f"- Additional measurements recorded as evidence: **{informational}**",
        f"- Scripts completing successfully: **{sum(item['passed'] for item in execution['executions'])}/{len(execution['executions'])}**",
        "",
        "## Software versions",
        "",
        "| Package | Version(s) |",
        "| --- | --- |",
        *[
            f"| {package} | {', '.join(sorted(package_versions))} |"
            for package, package_versions in sorted(versions.items())
        ],
        "",
        "![Pass rate by comparison group](method-pass-rate.png)",
        "",
        "## Method groups",
        "",
        "| Method group | Reference | Dataset | Tier | Result | Metrics |",
        "| --- | --- | --- | --- | --- | ---: |",
    ]
    for record in summary["records"]:
        gated = gating_metrics(record)
        record_passed = sum(metric["passed"] for metric in gated)
        lines.append(
            f"| `{record['method']}` | {record['reference_package']} | {record['dataset']} | {record['tier']} | "
            f"{status(record['passed'])} | {record_passed}/{len(gated)} |"
        )

    failed_metrics = [
        (record["method"], metric)
        for record in summary["records"]
        for metric in gating_metrics(record)
        if not metric["passed"]
    ]
    lines.extend(["", "## Failed metrics", ""])
    if failed_metrics:
        lines.extend(
            [
                "| Method group | Metric | Observed | Criterion |",
                "| --- | --- | ---: | --- |",
            ]
        )
        for method, metric in failed_metrics:
            lines.append(f"| `{method}` | `{metric['metric']}` | {metric['observed']:.8g} | {metric['criterion']} |")
    else:
        lines.append("None.")

    recorded = [
        (record["method"], metric)
        for record in summary["records"]
        for metric in record["metrics"]
        if not is_gating(metric)
    ]
    if recorded:
        lines.extend(
            [
                "",
                "## Recorded measurements (not gating)",
                "",
                "These quantify behaviour rather than test CPU/GPU equivalence, so they are"
                " reported without deciding the outcome.",
                "",
                "| Method group | Measurement | Observed | Why it is not a criterion |",
                "| --- | --- | ---: | --- |",
            ]
        )
        for method, metric in recorded:
            lines.append(
                f"| `{method}` | `{metric['metric']}` | {metric['observed']:.8g} | {metric.get('basis', '')} |"
            )

    justified = [
        (record["method"], metric)
        for record in summary["records"]
        for metric in gating_metrics(record)
        if metric.get("basis")
    ]
    if justified:
        lines.extend(
            [
                "",
                "## Basis for gating criteria",
                "",
                "| Method group | Criterion | Threshold | Basis |",
                "| --- | --- | --- | --- |",
            ]
        )
        for method, metric in justified:
            lines.append(f"| `{method}` | `{metric['metric']}` | {metric['criterion']} | {metric['basis']} |")

    incomplete = [
        item for item in execution["executions"] if not item["passed"] or not item.get("result_present", True)
    ]
    lines.extend(["", "## Incomplete or failing scripts", ""])
    if incomplete:
        lines.extend(["| Script | Exit code | Result record | Log |", "| --- | ---: | --- | --- |"])
        for item in incomplete:
            lines.append(
                f"| `{item['script']}` | {item['exit_code']} | {'yes' if item.get('result_present', True) else 'no'} | "
                f"[{item.get('log', 'not captured')}]({item.get('log', '#')}) |"
            )
    else:
        lines.append("None.")

    lines.extend(
        [
            "",
            "## Reviewer-facing evidence",
            "",
            "| Reviewer request | Evidence in this report |",
            "| --- | --- |",
            "| Numerical equivalence for deterministic methods | Absolute error and correlation for preprocessing, HVG, PCA, spatial statistics, activity inference, and perturbation methods |",
            "| Biological equivalence for stochastic methods | Label-invariant ARI/NMI, embedding trustworthiness, and neighborhood overlap |",
            "| Explicit clustering agreement | Leiden, Louvain, k-means, and spatial-niche ARI/NMI |",
            "| Explicit marker preservation | Per-cell-type top-50 overlap and score agreement |",
            "| Explicit embedding comparison | Quantitative overlap plus side-by-side biological-pipeline UMAPs |",
            "| Cell-type interpretation | Held-out annotation accuracy and CPU/GPU prediction agreement |",
            "| Additional scverse APIs | Direct Squidpy, Decoupler, and Pertpy reference comparisons |",
            "",
            "Raw metric records are available in [`metrics.csv`](metrics.csv) and [`equivalence.json`](../equivalence.json).",
            "Automatic publication of this report requires GPU-backed CI, either through a self-hosted GPU runner or a CI-to-Slurm integration.",
        ]
    )
    output.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(args.summary.read_text())
    execution = json.loads(args.execution.read_text())
    render_csv(summary, args.output_dir / "metrics.csv")
    render_plot(summary, args.output_dir / "method-pass-rate.png")
    render_markdown(summary, execution, args.output_dir / "summary.md")
    print(f"Wrote reviewer report to {args.output_dir}")


if __name__ == "__main__":
    main()
