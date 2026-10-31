# Legacy standalone scripts

These single-file scripts predate the `h2multisim` package and are kept for
reference and reproducibility of earlier results only. They are not imported,
tested, linted or maintained, and do not include the fixes made since v1.1
(see `CHANGELOG.md`).

Use the package instead:

```bash
pip install -e .
h2multisim run --scale 0d --scenario teesside
h2multisim run --config examples/04_config_driven_2d.json
```

`key_improvements_in_v1.3.txt` describes `h2_storage_sim (v1.3).py`, not the
`h2multisim` 1.3.x package.
