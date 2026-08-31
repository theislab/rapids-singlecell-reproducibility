from __future__ import annotations

import sys
from functools import partial
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))

import report
from report import capture, measure, package_version  # re-exported for this group's scripts

write_report = partial(report.write_report, reference_package="squidpy")
