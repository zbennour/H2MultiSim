# Security Policy

## Supported versions

Security fixes are made on the latest minor release of `h2multisim` on the `main` branch.

## Reporting a vulnerability

Please report suspected vulnerabilities privately through GitHub's
[private vulnerability reporting](https://github.com/zbennour/H2MultiSim/security/advisories/new)
rather than in a public issue. Include the affected version, a description of the issue and steps
to reproduce. You should receive an acknowledgement within a week.

H2MultiSim has no runtime dependencies and performs no network access. Its main input surface is
JSON/YAML run configuration files (YAML is parsed with `yaml.safe_load`) and the paths the CLI
writes outputs to.
