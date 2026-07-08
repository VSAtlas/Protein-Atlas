# Private Release Checklist (BYOL / License Safety)

This checklist is intended for internal release prep and is not part of the
public README surface.

## License and Redistribution

- Confirm root `LICENSE` is present and matches intended project licensing.
- Verify no proprietary binaries are committed to the repo.
- Verify no license keys, tokens, or license files are committed.
- Verify no private/internal download URLs are present in install scripts/docs.

## Installer and Environment Posture

- Confirm `tools/install_atlas_dev.sh` defaults to `environment.core.yml`.
- Confirm optional extras are install-on-demand (`ATLAS_INSTALL_OPTIONAL=1`).
- Confirm docs state BYOL for licensed/proprietary external tools.
- Confirm legacy `environment.yml` is clearly marked as compatibility-only.

## Tool Verification Expectations

- Run `atlas --verify-tools` in a configured environment.
- Confirm success criteria:
  - exit code `0`;
  - `Tools all successfully verified` appears;
  - runtime probe lines appear for Meeko ligand prep.
- Confirm optional missing tools do not fail verification unexpectedly.

## Public Smoke Path

- Run `bash tools/public_smoke_check.sh` in a fresh clone.
- Confirm script exits `0` without requiring proprietary tool installs.
- Confirm CI `cli-public-smoke` workflow path is green.

## Final Compliance Sweep

- Run secret scanning and dependency/license checks used by your org.
- Re-check `.gitignore` coverage for generated artifacts and local credentials.
- Capture a short release note with installation constraints and BYOL language.
