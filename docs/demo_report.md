# Demo Report Fixture

The deterministic report-demo fixture lives at:

```bash
chemdb/tests/fixtures/report_demo/heatmap_input.csv
```

It is synthetic and does not require docking tools. The fixture contains three ligands, three targets, one strong hit, weak/background rows, and an invalid pose marker.

Regenerate the test report in CI/local validation with:

```bash
pytest chemdb/tests/test_report_demo_fixture.py -q
```

The test writes generated HTML under pytest's temporary directory. Do not commit generated `report.html` unless it is intentionally added under a small `docs/examples/` artifact path.
