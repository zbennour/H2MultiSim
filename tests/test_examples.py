"""Every shipped example must run (scripts) or load-and-run (JSON configs) cleanly."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
ENV = {**os.environ, "PYTHONPATH": str(ROOT / "src")}


@pytest.mark.parametrize("script", sorted(p.name for p in EXAMPLES.glob("*.py")))
def test_example_script_runs(script: str, tmp_path: pathlib.Path) -> None:
    proc = subprocess.run(
        [sys.executable, str(EXAMPLES / script)],
        cwd=tmp_path,
        env=ENV,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]


@pytest.mark.parametrize("config", sorted(p.name for p in EXAMPLES.glob("*.json")))
def test_example_config_runs(config: str, tmp_path: pathlib.Path) -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "h2multisim", "run", "--config", str(EXAMPLES / config),
         "--out-dir", str(tmp_path), "--json", "-"],
        env=ENV, capture_output=True, text=True, timeout=300, check=False,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert json.loads(proc.stdout)["metrics"]["converged"] == 1.0
