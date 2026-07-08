# Development

## Pre-commit

The repository includes a conservative pre-commit configuration for basic file
hygiene and Ruff checks on maintained code paths. It intentionally does not run
the full test suite or MyPy on every commit.

Install and run locally:

```bash
pre-commit install
pre-commit run --all-files
```

The hooks check trailing whitespace, final newlines, YAML/TOML syntax, merge
conflict markers, and Ruff autofixes for maintained Python code.

## Before Opening a PR

Run the focused tests for the area you touched, then run:

```bash
atlas dev verify --full --fix --smoke
```

When environment activation is uncertain, run:

```bash
atlas --doctor
```
