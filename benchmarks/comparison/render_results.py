from __future__ import annotations

import csv
import json
from datetime import UTC, datetime
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path


def status(value: bool) -> str:
    return "PASS" if value else "FAIL"


def cuda_version(packed: int) -> str:
    """12090 -> "12.9"."""
    return f"{packed // 1000}.{packed % 1000 // 10}"


def gpu_provenance(hardware: Path) -> list[str]:
    """Hardware rows for the report."""
    if not hardware.exists():
        return []
    recorded = json.loads(hardware.read_text())
    device = recorded["device"]
    return [
        "## Hardware",
        "",
        "| | |",
        "| --- | --- |",
        f"| GPU | {device['device_name']} |",
        f"| Compute capability | {device['compute_capability']} |",
        f"| CUDA driver / runtime | {cuda_version(device['driver_version'])} / "
        f"{cuda_version(device['cuda_runtime_version'])} |",
        f"| Device memory | {device['total_memory_bytes'] / 1024**3:.1f} GiB |",
        f"| Measured | {recorded['measured_at']} |",
        "",
    ]


def scale(record: dict) -> str:
    """The compared input as `cells x features`, from the record the run wrote."""
    shape = record.get("shape")
    return " x ".join(f"{dimension:,}" for dimension in shape) if shape else ""


def gating_metrics(record: dict) -> list[dict]:
    return [metric for metric in record["metrics"] if metric["gating"]]


def render_csv(summary: dict, output: Path) -> None:
    rows = [
        {
            "method": record["method"],
            "reference_package": record["reference_package"],
            "dataset": record["dataset"],
            "shape": scale(record),
            "metric": metric["metric"],
            "observed": metric["observed"],
            "criterion": metric["criterion"],
            "passed": metric["passed"],
            "versions": json.dumps(record["versions"], sort_keys=True),
        }
        for record in summary["records"]
        for metric in record["metrics"]
    ]
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def render_markdown(summary: dict, execution: dict, output: Path, hardware: Path) -> None:
    versions = {}
    for record in summary["records"]:
        for package, package_version in record["versions"].items():
            versions.setdefault(package, set()).add(str(package_version))
    lines = [
        "# CPU/GPU equivalence report",
        "",
        f"Scored {datetime.now(UTC):%Y-%m-%d}.",
        "",
        "## Outcome",
        "",
        *(
            [
                "> Partial run: records for comparisons not rerun come from an earlier run.",
                "",
            ]
            if execution.get("partial")
            else []
        ),
        f"- Overall: **{status(summary['passed'])}**",
        f"- Gating metrics passing: **{summary['n_passed_metrics']}/{summary['n_metrics']}**",
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
        *gpu_provenance(hardware),
        "## Method groups",
        "",
        "| Method group | Reference | Dataset | Scale | Result | Metrics |",
        "| --- | --- | --- | --- | --- | ---: |",
    ]
    for record in summary["records"]:
        gated = gating_metrics(record)
        record_passed = sum(metric["passed"] for metric in gated)
        lines.append(
            f"| `{record['method']}` | {record['reference_package']} | {record['dataset']} | "
            f"{scale(record)} | {status(record['passed'])} | {record_passed}/{len(gated)} |"
        )

    gated = sorted(
        ((record["method"], metric) for record in summary["records"] for metric in gating_metrics(record)),
        key=lambda item: item[1]["passed"],
    )
    lines.extend(
        [
            "",
            "## Criteria",
            "",
            "| Method group | Metric | Observed | Criterion | Result |",
            "| --- | --- | ---: | --- | --- |",
            *[
                f"| `{method}` | `{metric['metric']}` | {metric['observed']:.8g} | {metric['criterion']} | "
                f"{status(metric['passed'])} |"
                for method, metric in gated
            ],
        ]
    )

    recorded = [
        (record["method"], metric)
        for record in summary["records"]
        for metric in record["metrics"]
        if not metric["gating"]
    ]
    if recorded:
        lines.extend(
            [
                "",
                "## Recorded measurements (not gating)",
                "",
                "| Method group | Measurement | Observed |",
                "| --- | --- | ---: |",
            ]
        )
        for method, metric in recorded:
            lines.append(f"| `{method}` | `{metric['metric']}` | {metric['observed']:.8g} |")

    incomplete = [item for item in execution["executions"] if not item["passed"]]
    if incomplete:
        lines.extend(["", "## Incomplete or failing scripts", ""])
        lines.extend(["| Script | Exit code | Result record | Log |", "| --- | ---: | --- | --- |"])
        for item in incomplete:
            lines.append(
                f"| `{item['script']}` | {item['exit_code']} | {'yes' if item['result_present'] else 'no'} | "
                f"[{item['log']}]({item['log']}) |"
            )

    output.write_text("\n".join(lines) + "\n")


def render(summary: dict, execution: dict, report_dir: Path, hardware: Path) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)
    render_csv(summary, report_dir / "metrics.csv")
    render_markdown(summary, execution, report_dir / "summary.md", hardware)
