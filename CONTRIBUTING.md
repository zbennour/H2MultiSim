# Contributing to H2MultiSim

Thank you for considering a contribution. H2MultiSim is a small, dependency-free
scientific Python project; please keep changes focused and tested.

## Getting started

```bash
git clone https://github.com/zbennour/H2MultiSim.git
cd H2MultiSim
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Development workflow

- Create a branch from `main`: `git checkout -b feat/my-change`.
- Run the full suite before pushing: `pytest` (32 tests, ~65 s).
- Keep `ruff check .` and `mypy src` clean.
- Add or update tests for any behaviour change; prefer tests that close a mass
  balance or compare against a reference (CoolProp, analytic, or an independent
  simulator) rather than sign/bound checks.
- Update `CHANGELOG.md` and the README benchmark table if you change physics.

## Scope and style

- Runtime dependencies are not allowed in the core package. Optional extras
  (e.g. PyYAML) must be guarded with `try/except ImportError`.
- Keep the public API stable; new options should be keyword-only with defaults.
- Physics changes need a short derivation or reference in the PR description.

## Submitting

Open a pull request against `main`. CI runs the test suite on Python 3.9–3.13.
A maintainer will review; expect discussion on numerical accuracy and mass balance.

## Reporting issues

Please include: Python version, H2MultiSim version (`h2multisim --version`),
the config or scenario, and the full traceback. For physics discrepancies, a
CoolProp/REFPROP reference value is invaluable.
