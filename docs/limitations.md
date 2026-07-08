## Atlas2 Limitations

- Atlas outputs are computational hypotheses, not validated binding claims.
- Docking scores are not experimental affinity measurements.
- Atlas score is a relative plausibility/ranking signal against the selected reference/diversity library.
- Atlas score, `z_selected`, and related reporting fields should be interpreted relative to the configured comparator/reference library, not as universal binding constants.
- Protein preparation choices can materially change outcomes, including protonation state, missing loops, alternate locations, chain selection, retained waters, ions/metals, binding-site definition, and APO/HOLO structure choice.
- Ligand preparation choices can materially change outcomes, including tautomer/protomer state, stereochemistry, salt handling, conformer generation, and pH-specific enumeration.
- MM/GBSA rescoring is supportive evidence, not definitive proof of binding.
- ADR linkage fields require independent biological and literature validation.
- Reported rankings are sensitive to input library composition and target structure quality.
- Side-effect and ADR hypotheses require independent evidence for exposure, tissue expression, pathway relevance, target engagement, and clinical/biological plausibility.
- A favorable docking or rescoring result does not prove causality for an adverse event.
- Benchmark or DUD/decoy performance does not guarantee prospective performance on new targets or chemotypes.
- External tool versions, hardware, thread scheduling, and distributed execution can introduce small reproducibility differences.
- Atlas reports should state the input structures, ligand libraries, tool versions, and scoring/rescoring paths used for any publication claim.
