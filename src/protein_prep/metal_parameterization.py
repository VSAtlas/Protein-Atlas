"""Metal force-field parameterization planning helpers."""

from __future__ import annotations

import json
import os
import shutil
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence, cast

from protein_prep.pdb_fixer_runtime import (
    load_canonical_cofactors,
    load_canonical_metals,
    load_canonical_waters,
)
from protein_prep.geometry.repair import repair_terminal_oxt_geometry_in_pdb


MCPB_TOOLS: tuple[str, ...] = (
    "MCPB.py",
    "tleap",
    "antechamber",
    "parmchk2",
    "metalpdb2mol2.py",
    "pdb4amber",
    "cpptraj",
    "parmed",
    "rungms",
    "g16",
    "g09",
    "g03",
    "formchk",
    "sqm",
)

STANDARD_RESNAMES: set[str] = {
    "ALA",
    "ARG",
    "ASH",
    "ASN",
    "ASP",
    "CYM",
    "CYS",
    "GLH",
    "GLN",
    "GLU",
    "GLY",
    "HID",
    "HIE",
    "HIP",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
WATER_RESNAMES: set[str] = {"HOH", "WAT", "DOD", "TIP", "SOL"}
METAL_RESNAMES: set[str] = {
    "AG",
    "AL",
    "BA",
    "CA",
    "CD",
    "CO",
    "CR",
    "CS",
    "CU",
    "FE",
    "HG",
    "K",
    "LI",
    "MG",
    "MN",
    "MO",
    "NA",
    "NI",
    "RB",
    "SR",
    "V",
    "W",
    "YB",
    "ZN",
}
MCPB_AUTOPARAM_ENV = "ATLAS_MCPB_AUTOPARAM_COFACTORS"


def build_mcpb_parameterization_plan(
    *,
    receptor_pdb: Path | None,
    metal_rows: Sequence[Mapping[str, object]],
    work_dir: Path,
) -> dict[str, object]:
    """Build an auditable AmberTools MCPB.py handoff plan for metal centers."""

    sites = _parameterized_sites(metal_rows)
    tools = _tool_paths(MCPB_TOOLS)
    missing_required = _missing_tools(tools, ("MCPB.py", "tleap"))
    amberhome = _amberhome_from_tools(tools)
    work_dir.mkdir(parents=True, exist_ok=True)
    alias_sets = _mcpb_alias_sets()
    cofactor_plan = _mcpb_cofactor_parameterization_plan(
        receptor_pdb=receptor_pdb,
        sites=sites,
        work_dir=work_dir,
        tools=tools,
        amberhome=amberhome,
        alias_sets=alias_sets,
    )
    preclean_plan = _mcpb_preclean_plan(sites, alias_sets, cofactor_plan)
    receptor_info = _write_mcpb_receptor(work_dir, receptor_pdb, preclean_plan)
    templates = _write_mcpb_templates(work_dir, sites, receptor_info, cofactor_plan)
    status = _parameterization_status(
        site_count=len(sites),
        missing_required=missing_required,
        receptor_ready=receptor_pdb is not None and receptor_pdb.exists(),
        templates=templates,
    )
    return {
        "status": status,
        "site_count": len(sites),
        "missing_required_tools": missing_required,
        "optional_missing_tools": _optional_missing_tools(tools, missing_required),
        "tool_paths": tools,
        "amberhome": str(amberhome) if amberhome else "",
        "receptor_pdb": str(receptor_pdb) if receptor_pdb else "",
        "mcpb_receptor_pdb": str(receptor_info.get("path", "")),
        "mcpb_preclean": receptor_info.get("preclean", {}),
        "mcpb_cofactor_parameterization": cofactor_plan,
        "work_dir": str(work_dir),
        "treatment_counts": _format_counter(
            Counter(str(site.get("metal_treatment", "")) for site in sites)
        ),
        "templates": templates,
        "commands": _mcpb_commands(templates, tools, amberhome),
        "note": (
            "MCPB.py requires curated metal-site residue selections, residue/atom naming, "
            "metal oxidation state, ligand/cofactor mol2/frcmod files, and QM/RESP or "
            "empirical follow-through. Atlas pre-cleans only the MCPB staging copy by "
            "removing hydrogens, normalizing MCPB residue identifiers, and applying "
            "unambiguous local metal-bound residue-name heuristics. Canonical "
            "metal-bound cofactors are parameterized with the MMGBSA AmberTools "
            "fallback only when an SDF is already available; this is not an automatic "
            "force-field guarantee."
        ),
    }


def summarize_parameterization_plan(plan: Mapping[str, object]) -> dict[str, object]:
    """Flatten a parameterization plan for benchmark CSV output."""

    missing = plan.get("missing_required_tools", [])
    optional_missing = plan.get("optional_missing_tools", [])
    return {
        "metal_parameterization_status": str(plan.get("status", "")),
        "metal_parameterization_site_count": _safe_int(plan.get("site_count", 0)),
        "metal_parameterization_missing_tools": ",".join(missing)
        if isinstance(missing, list)
        else str(missing or ""),
        "metal_parameterization_optional_missing_tools": ",".join(optional_missing)
        if isinstance(optional_missing, list)
        else str(optional_missing or ""),
        **_summarize_preclean(plan.get("mcpb_preclean", {})),
        **_summarize_cofactor_plan(plan.get("mcpb_cofactor_parameterization", {})),
    }


def write_parameterization_plan(path: Path, plan: Mapping[str, object]) -> None:
    """Write a parameterization sidecar."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(plan, indent=2, sort_keys=True), encoding="utf-8")


def _tool_paths(names: Sequence[str]) -> dict[str, str]:
    return {name: _resolve_ambertools_tool(name) for name in names}


def _resolve_ambertools_tool(name: str) -> str:
    found = shutil.which(name)
    if found:
        return found
    for prefix in _ambertools_prefix_candidates():
        candidate = prefix / "bin" / name
        if candidate.is_file():
            return str(candidate)
    return ""


def _ambertools_prefix_candidates() -> list[Path]:
    candidates: list[Path] = []
    for raw in (
        os.environ.get("MMGBSA_AMBERTOOLS_PREFIX"),
        os.environ.get("AMBERTOOLS_PREFIX"),
        os.environ.get("AMBERHOME"),
    ):
        if raw:
            candidates.append(Path(raw).expanduser())
    repo_path = Path(__file__).resolve()
    for parent in repo_path.parents:
        candidates.append(parent / "tools" / "envs" / "ambertools")
    home = Path.home()
    candidates.extend(
        [
            home / "tools" / "envs" / "ambertools",
            home / "micromamba" / "envs" / "AmberTools25",
        ]
    )
    return _dedup_existing_prefixes(candidates)


def _dedup_existing_prefixes(candidates: Sequence[Path]) -> list[Path]:
    seen: set[str] = set()
    prefixes: list[Path] = []
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        prefixes.append(candidate)
    return prefixes


def _amberhome_from_tools(tools: Mapping[str, str]) -> Path | None:
    mcpb = tools.get("MCPB.py", "")
    if not mcpb:
        return None
    path = Path(mcpb)
    if path.parent.name == "bin":
        return path.parent.parent
    return None


def _parameterized_sites(
    metal_rows: Sequence[Mapping[str, object]],
) -> list[Mapping[str, object]]:
    return [
        row
        for row in metal_rows
        if bool(row.get("requires_parameterization"))
        and row.get("metal_treatment") != "autodock4zn_candidate"
    ]


def _missing_tools(tools: Mapping[str, str], required: Sequence[str]) -> list[str]:
    return [name for name in required if not tools.get(name)]


def _optional_missing_tools(
    tools: Mapping[str, str],
    required_missing: Sequence[str],
) -> list[str]:
    required = set(required_missing)
    return [name for name, path in tools.items() if not path and name not in required]


def _parameterization_status(
    *,
    site_count: int,
    missing_required: Sequence[str],
    receptor_ready: bool,
    templates: Sequence[Mapping[str, object]],
) -> str:
    if site_count == 0:
        return "not_applicable"
    if not receptor_ready:
        return "missing_receptor"
    if missing_required:
        return "missing_tools"
    if any(bool(template.get("needs_manual_ids")) for template in templates):
        return "ready_for_curated_mcpb_inputs"
    return "ready_for_heuristic_mcpb_trial"


def _write_mcpb_templates(
    work_dir: Path,
    sites: Sequence[Mapping[str, object]],
    receptor_info: Mapping[str, object],
    cofactor_plan: Mapping[str, object],
) -> list[dict[str, object]]:
    work_dir.mkdir(parents=True, exist_ok=True)
    return [
        _write_mcpb_template(work_dir, idx, site, receptor_info, cofactor_plan)
        for idx, site in enumerate(sites, start=1)
    ]


def _write_mcpb_template(
    work_dir: Path,
    index: int,
    site: Mapping[str, object],
    receptor_info: Mapping[str, object],
    cofactor_plan: Mapping[str, object],
) -> dict[str, object]:
    site_id = _safe_site_id(site, index)
    template_path = work_dir / f"{site_id}.mcpb.in"
    metal_charge = _mcpb_charge(site)
    donors = site.get("donors", [])
    donor_ids = _donor_ids(donors if isinstance(donors, list) else [])
    protonation_notes = _protonation_notes(site)
    preclean_notes = _preclean_notes(receptor_info, site)
    ion_id = _ion_id(site, receptor_info)
    needs_manual_ids = ion_id == "REVIEW_METAL_ATOM_ID"
    metal_pdb = _write_metal_pdb(work_dir, site_id, site)
    metal_mol2 = work_dir / f"{site_id}.mol2"
    cofactor_entries = _stage_mcpb_cofactor_files(
        work_dir,
        _site_cofactor_parameter_entries(site, cofactor_plan),
    )
    naa_mol2files = _cofactor_files(cofactor_entries, "mcpb_mol2_token")
    frcmod_files = _cofactor_files(cofactor_entries, "mcpb_frcmod_token")
    lines = [
        "# Atlas-generated MCPB.py heuristic input. Review before publication use.",
        "# Confirm ion_ids, residue names, ligand/cofactor mol2/frcmod files, charge, spin,",
        "# multiplicity, QM/RESP settings, and final topology quality.",
        f"# Atlas metal site: {site_id}",
        f"# Element: {site.get('element', '')}",
        f"# Formal charge guess: {metal_charge}",
        f"# Coordination geometry: {site.get('coordination_geometry', '')}",
        f"# First-shell donors: {','.join(donor_ids) or 'none'}",
        f"# Protonation review: {'; '.join(protonation_notes) or 'none'}",
        f"# MCPB pre-clean: {'; '.join(preclean_notes) or 'none'}",
        f"original_pdb {str(receptor_info.get('path') or 'REVIEW_RECEPTOR.pdb')}",
        f"group_name {site_id}",
        f"ion_ids {ion_id}",
        f"ion_info {_mcpb_ion_info(site, metal_charge)}",
        f"ion_mol2files {metal_mol2.name}",
        _mcpb_file_line("naa_mol2files", naa_mol2files),
        _mcpb_file_line("frcmod_files", frcmod_files),
        "cut_off 2.8",
        "large_opt 1",
        "software_version gms",
    ]
    template_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "site_id": site_id,
        "template": str(template_path),
        "element": str(site.get("element", "")),
        "formal_charge_primary": str(site.get("formal_charge_primary", "")),
        "coordination_geometry": str(site.get("coordination_geometry", "")),
        "donor_ids": donor_ids,
        "ion_id": ion_id,
        "needs_manual_ids": needs_manual_ids,
        "metal_pdb": str(metal_pdb) if metal_pdb else "",
        "metal_mol2": str(metal_mol2),
        "metal_charge": metal_charge,
        "ion_info": _mcpb_ion_info(site, metal_charge),
        "protonation_notes": protonation_notes,
        "preclean_notes": preclean_notes,
        "cofactor_parameters": cofactor_entries,
        "naa_mol2files": naa_mol2files,
        "frcmod_files": frcmod_files,
        "naa_mol2file_paths": _cofactor_files(cofactor_entries, "mcpb_mol2_path"),
        "frcmod_file_paths": _cofactor_files(cofactor_entries, "mcpb_frcmod_path"),
        "stage1_expected_outputs": _stage1_expected_outputs(work_dir, site_id),
        "qm_required_outputs": _qm_required_outputs(work_dir, site_id),
    }


def _mcpb_commands(
    templates: Sequence[Mapping[str, object]],
    tools: Mapping[str, str],
    amberhome: Path | None,
) -> list[dict[str, object]]:
    commands: list[dict[str, object]] = []
    for template in templates:
        input_file = str(template.get("template", ""))
        if not input_file:
            continue
        metal_pdb = str(template.get("metal_pdb", ""))
        metal_mol2 = str(template.get("metal_mol2", ""))
        if metal_pdb and metal_mol2 and tools.get("metalpdb2mol2.py"):
            commands.append(
                {
                    "name": "metal_mol2_from_pdb",
                    "argv": _amber_argv(
                        tools,
                        amberhome,
                        [
                            tools["metalpdb2mol2.py"],
                            "-i",
                            metal_pdb,
                            "-o",
                            metal_mol2,
                            "-c",
                            str(template.get("metal_charge", "REVIEW")),
                        ],
                    ),
                }
            )
        commands.extend(
            [
                {
                    "name": "mcpb_stage_1_generate_models",
                    "argv": _mcpb_argv(tools, amberhome, input_file, "1"),
                },
                *_mcpb_qm_commands(template, tools),
                {
                    "name": "mcpb_stage_2_parameterize",
                    "argv": _mcpb_argv(tools, amberhome, input_file, "2"),
                    "requires_existing_files": _stage_2_qm_requirements(template),
                },
                {
                    "name": "mcpb_stage_3_resp",
                    "argv": _mcpb_argv(tools, amberhome, input_file, "3"),
                    "requires_existing_files": _stage_3_qm_requirements(template),
                },
                {
                    "name": "mcpb_stage_4_tleap",
                    "argv": _mcpb_argv(tools, amberhome, input_file, "4"),
                },
                {
                    "name": "mcpb_stage_4_nonbonded_12_6",
                    "argv": _mcpb_argv(tools, amberhome, input_file, "4n2"),
                    "note": "Nonbonded model fallback; requires ion_info and is not a bonded MCPB replacement.",
                },
            ]
        )
    return commands


def _mcpb_argv(
    tools: Mapping[str, str],
    amberhome: Path | None,
    input_file: str,
    step: str,
) -> list[str]:
    return _amber_argv(
        tools,
        amberhome,
        [tools.get("MCPB.py", "MCPB.py"), "-i", input_file, "-s", step],
    )


def _amber_argv(
    tools: Mapping[str, str],
    amberhome: Path | None,
    argv: Sequence[str],
) -> list[str]:
    if amberhome:
        return ["env", f"AMBERHOME={amberhome}", *argv]
    return list(argv)


def _mcpb_file_line(key: str, paths: Sequence[str]) -> str:
    if paths:
        return f"{key} {' '.join(paths)}"
    suffix = "mol2" if key == "naa_mol2files" else "frcmod"
    return f"# {key} REVIEW_NONSTANDARD_RESIDUES.{suffix}"


def _mcpb_qm_commands(
    template: Mapping[str, object],
    tools: Mapping[str, str],
) -> list[dict[str, object]]:
    work_dir = Path(str(template.get("template", ""))).parent
    site_id = str(template.get("site_id", ""))
    if not site_id:
        return []
    commands: list[dict[str, object]] = []
    if tools.get("rungms"):
        commands.extend(_gamess_qm_commands(work_dir, site_id, tools["rungms"]))
    gaussian = _gaussian_tool(tools)
    if gaussian:
        commands.extend(_gaussian_qm_commands(work_dir, site_id, gaussian, tools))
    if not commands:
        commands.append(
            {
                "name": "mcpb_qm_required_manual",
                "argv": [],
                "status": "missing_qm_engine",
                "note": (
                    "MCPB stage 2/3 require Gaussian or GAMESS-US output logs. "
                    "GAMESS-US is free but not open-source/redistributable; install "
                    "or module-load rungms externally, or provide equivalent logs."
                ),
            }
        )
    return commands


def _gamess_qm_commands(
    work_dir: Path,
    site_id: str,
    rungms: str,
) -> list[dict[str, object]]:
    return [
        {
            "name": "mcpb_qm_gamess_small_opt",
            "argv": [rungms, f"{site_id}_small_opt", "01", "2", "1"],
            "cwd": str(work_dir),
            "stdout": str(work_dir / f"{site_id}_small_opt.log"),
            "note": (
                "After this finishes, copy optimized coordinates from the log into "
                f"{site_id}_small_fc.inp before running the force-constant job."
            ),
        },
        {
            "name": "mcpb_qm_gamess_small_fc",
            "argv": [rungms, f"{site_id}_small_fc", "01", "2", "1"],
            "cwd": str(work_dir),
            "stdout": str(work_dir / f"{site_id}_small_fc.log"),
        },
        {
            "name": "mcpb_qm_gamess_large_mk",
            "argv": [rungms, f"{site_id}_large_mk", "01", "2", "1"],
            "cwd": str(work_dir),
            "stdout": str(work_dir / f"{site_id}_large_mk.log"),
        },
    ]


def _gaussian_qm_commands(
    work_dir: Path,
    site_id: str,
    gaussian: str,
    tools: Mapping[str, str],
) -> list[dict[str, object]]:
    commands = [
        _stdin_stdout_command(
            "mcpb_qm_gaussian_small_opt",
            gaussian,
            work_dir / f"{site_id}_small_opt.com",
            work_dir / f"{site_id}_small_opt.log",
        ),
        _stdin_stdout_command(
            "mcpb_qm_gaussian_small_fc",
            gaussian,
            work_dir / f"{site_id}_small_fc.com",
            work_dir / f"{site_id}_small_fc.log",
        ),
        _stdin_stdout_command(
            "mcpb_qm_gaussian_large_mk",
            gaussian,
            work_dir / f"{site_id}_large_mk.com",
            work_dir / f"{site_id}_large_mk.log",
        ),
    ]
    if tools.get("formchk"):
        commands.append(
            {
                "name": "mcpb_qm_gaussian_formchk",
                "argv": [
                    tools["formchk"],
                    str(work_dir / f"{site_id}_small_opt.chk"),
                    str(work_dir / f"{site_id}_small_opt.fchk"),
                ],
                "cwd": str(work_dir),
            }
        )
    return commands


def _stdin_stdout_command(name: str, exe: str, stdin: Path, stdout: Path) -> dict[str, object]:
    return {
        "name": name,
        "argv": [exe],
        "stdin": str(stdin),
        "stdout": str(stdout),
        "cwd": str(stdin.parent),
    }


def _gaussian_tool(tools: Mapping[str, str]) -> str:
    for name in ("g16", "g09", "g03"):
        if tools.get(name):
            return tools[name]
    return ""


def _stage1_expected_outputs(work_dir: Path, site_id: str) -> dict[str, str]:
    return {
        "small_pdb": str(work_dir / f"{site_id}_small.pdb"),
        "standard_pdb": str(work_dir / f"{site_id}_standard.pdb"),
        "large_pdb": str(work_dir / f"{site_id}_large.pdb"),
        "small_opt_input": str(work_dir / f"{site_id}_small_opt.inp"),
        "small_fc_input": str(work_dir / f"{site_id}_small_fc.inp"),
        "large_mk_input": str(work_dir / f"{site_id}_large_mk.inp"),
    }


def _qm_required_outputs(work_dir: Path, site_id: str) -> dict[str, str]:
    return {
        "small_force_constant_log": str(work_dir / f"{site_id}_small_fc.log"),
        "large_mk_esp_log": str(work_dir / f"{site_id}_large_mk.log"),
        "small_opt_formchk_gaussian_only": str(work_dir / f"{site_id}_small_opt.fchk"),
    }


def _stage_2_qm_requirements(template: Mapping[str, object]) -> list[str]:
    outputs = template.get("qm_required_outputs", {})
    if not isinstance(outputs, Mapping):
        return []
    required = [str(outputs.get("small_force_constant_log", ""))]
    gaussian_fchk = str(outputs.get("small_opt_formchk_gaussian_only", ""))
    if str(template.get("software_version", "gms")) != "gms" and gaussian_fchk:
        required.append(gaussian_fchk)
    return [path for path in required if path]


def _stage_3_qm_requirements(template: Mapping[str, object]) -> list[str]:
    outputs = template.get("qm_required_outputs", {})
    if not isinstance(outputs, Mapping):
        return []
    path = str(outputs.get("large_mk_esp_log", ""))
    return [path] if path else []


def _mcpb_alias_sets() -> dict[str, set[str]]:
    return {
        "metals": set(load_canonical_metals(None)) or set(METAL_RESNAMES),
        "waters": set(load_canonical_waters(None)) or set(WATER_RESNAMES),
        "cofactors": set(load_canonical_cofactors(None)),
    }


def _mcpb_cofactor_parameterization_plan(
    *,
    receptor_pdb: Path | None,
    sites: Sequence[Mapping[str, object]],
    work_dir: Path,
    tools: Mapping[str, str],
    amberhome: Path | None,
    alias_sets: Mapping[str, set[str]],
) -> dict[str, object]:
    residues = _metal_bound_cofactor_residues(sites, alias_sets)
    entries: list[dict[str, object]] = []
    if not residues:
        return _cofactor_plan_payload("not_applicable", entries, residues)

    for residue in residues:
        entries.append(
            _parameterize_cofactor_residue(
                receptor_pdb=receptor_pdb,
                residue=residue,
                work_dir=work_dir,
                tools=tools,
                amberhome=amberhome,
            )
        )
    status = "parameterized" if all(bool(e.get("ok")) for e in entries) else "partial"
    if not any(bool(e.get("ok")) for e in entries):
        status = "needs_curated_parameters"
    return _cofactor_plan_payload(status, entries, residues)


def _cofactor_plan_payload(
    status: str,
    entries: Sequence[Mapping[str, object]],
    residues: Sequence[Mapping[str, str]],
) -> dict[str, object]:
    return {
        "status": status,
        "autoparam_enabled": _mcpb_autoparam_enabled(),
        "cofactor_residue_count": len(residues),
        "parameterized_count": sum(1 for entry in entries if bool(entry.get("ok"))),
        "unresolved_count": sum(1 for entry in entries if not bool(entry.get("ok"))),
        "entries": list(entries),
        "note": (
            "Canonical metal-bound cofactors are parameterized with the existing "
            "MMGBSA AmberTools antechamber/parmchk2 fallback only when an SDF is "
            "already available near the benchmark/receptor directory. Missing or "
            "failed cofactors remain curated-parameter requirements."
        ),
    }


def _metal_bound_cofactor_residues(
    sites: Sequence[Mapping[str, object]],
    alias_sets: Mapping[str, set[str]],
) -> list[dict[str, str]]:
    cofactors = alias_sets.get("cofactors", set())
    residues: dict[tuple[str, str, str], dict[str, str]] = {}
    for site in sites:
        donors = site.get("donors", [])
        if not isinstance(donors, list):
            continue
        for donor in donors:
            if not isinstance(donor, Mapping):
                continue
            resname = str(donor.get("resname", "")).upper()
            if resname not in cofactors:
                continue
            key = (resname, str(donor.get("chain", "")), str(donor.get("resseq", "")))
            residues[key] = {
                "resname": key[0],
                "chain": key[1],
                "resseq": key[2],
            }
    return [residues[key] for key in sorted(residues)]


def _parameterize_cofactor_residue(
    *,
    receptor_pdb: Path | None,
    residue: Mapping[str, str],
    work_dir: Path,
    tools: Mapping[str, str],
    amberhome: Path | None,
) -> dict[str, object]:
    resname = residue["resname"]
    sdf = _find_cofactor_sdf(receptor_pdb, work_dir, resname)
    out_root = work_dir / "cofactor_parameters"
    mol2 = out_root / "mol2" / "mcpb" / f"{resname}.mol2"
    frcmod = out_root / "frcmod" / "mcpb" / f"{resname}.frcmod"
    mol2.parent.mkdir(parents=True, exist_ok=True)
    frcmod.parent.mkdir(parents=True, exist_ok=True)
    base = {
        **dict(residue),
        "sdf_path": str(sdf) if sdf else "",
        "mol2_path": str(mol2),
        "frcmod_path": str(frcmod),
    }
    if not sdf:
        return {**base, "ok": False, "status": "missing_sdf"}
    if not _mcpb_autoparam_enabled():
        return {**base, "ok": False, "status": "autoparam_disabled"}
    if not tools.get("antechamber") or not tools.get("parmchk2"):
        return {**base, "ok": False, "status": "missing_ambertools"}
    try:
        from post_docking.mmgbsa.prep_for_mmgbsa import (
            parameterize_ligand_with_fallback,
        )

        cfg = _mmgbsa_cofactor_cfg(amberhome)
        result = parameterize_ligand_with_fallback(
            sdf_path=str(sdf),
            out_mol2=str(mol2),
            out_frcmod=str(frcmod),
            cfg=cfg,
            stage_dir="mcpb",
            ligand_stem=resname,
            force=False,
        )
    except Exception as exc:
        return {**base, "ok": False, "status": "parameterization_exception", "error": str(exc)}
    atom_name_map = _cofactor_mol2_atom_name_map(receptor_pdb, residue, sdf, mol2)
    return {
        **base,
        "ok": bool(result.get("ok")),
        "status": "parameterized" if bool(result.get("ok")) else "parameterization_failed",
        "publication_ready": bool(result.get("publication_ready")),
        "charge_method": str(result.get("charge_method", "")),
        "net_charge_used": result.get("net_charge_used"),
        "metadata_path": str(result.get("metadata_path", "")),
        **atom_name_map,
    }


def _mcpb_autoparam_enabled() -> bool:
    raw = os.environ.get(MCPB_AUTOPARAM_ENV, "1").strip().lower()
    return raw not in {"0", "false", "no", "off"}


def _mmgbsa_cofactor_cfg(amberhome: Path | None) -> dict[str, object]:
    return {
        "MMGBSA_AMBERTOOLS_PREFIX": str(amberhome) if amberhome else "",
        "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE": False,
        "MMGBSA_INPUT_CHEMISTRY_SOURCE": "ccd_or_local_sdf_for_mcpb_cofactor",
        "MMGBSA_LIGAND_PRIMARY_CHARGE_METHOD": "bcc",
        "MMGBSA_LIGAND_FALLBACK_CHARGE_METHOD": "gas",
        "MMGBSA_LIGAND_AT": "gaff2",
        "MMGBSA_LIGAND_CHEMISTRY_STRICT": False,
        "MMGBSA_RDKit_VALIDATE": True,
    }


def _find_cofactor_sdf(
    receptor_pdb: Path | None,
    work_dir: Path,
    resname: str,
) -> Path | None:
    roots = _cofactor_sdf_roots(receptor_pdb, work_dir)
    names = (
        f"{resname}_ideal.sdf",
        f"{resname}.sdf",
        f"{resname.lower()}_ideal.sdf",
        f"{resname.lower()}.sdf",
    )
    for root in roots:
        for name in names:
            candidate = root / name
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
    return None


def _cofactor_sdf_roots(receptor_pdb: Path | None, work_dir: Path) -> list[Path]:
    candidates = [work_dir, work_dir.parent]
    if receptor_pdb is not None:
        candidates.extend([receptor_pdb.parent, receptor_pdb.parent.parent])
    seen: set[str] = set()
    roots: list[Path] = []
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        roots.append(candidate)
    return roots


def _site_cofactor_parameter_entries(
    site: Mapping[str, object],
    cofactor_plan: Mapping[str, object],
) -> list[dict[str, object]]:
    site_resnames = _site_nonprotein_resnames(site)
    entries = cofactor_plan.get("entries", [])
    if not isinstance(entries, list):
        return []
    return [
        dict(entry)
        for entry in entries
        if isinstance(entry, Mapping)
        and bool(entry.get("ok"))
        and str(entry.get("resname", "")).upper() in site_resnames
    ]


def _stage_mcpb_cofactor_files(
    work_dir: Path,
    entries: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    staged: list[dict[str, object]] = []
    for entry in entries:
        staged_entry = dict(entry)
        resname = str(entry.get("resname", "")).upper()
        mol2_target = work_dir / f"{resname}.mol2"
        frcmod_target = work_dir / f"{resname}.frcmod"
        try:
            _copy_mcpb_parameter_file(
                entry.get("mol2_path"),
                mol2_target,
                resname,
                entry.get("mol2_atom_name_map"),
            )
            _copy_mcpb_parameter_file(entry.get("frcmod_path"), frcmod_target, "", None)
        except OSError as exc:
            staged_entry.update(
                {
                    "mcpb_stage_ok": False,
                    "mcpb_stage_error": str(exc),
                }
            )
            staged.append(staged_entry)
            continue
        staged_entry.update(
            {
                "mcpb_stage_ok": True,
                "mcpb_mol2_path": str(mol2_target),
                "mcpb_frcmod_path": str(frcmod_target),
                "mcpb_mol2_token": mol2_target.name,
                "mcpb_frcmod_token": frcmod_target.name,
            }
        )
        staged.append(staged_entry)
    return staged


def _copy_mcpb_parameter_file(
    source: object,
    target: Path,
    residue_name: str,
    atom_name_map: object,
) -> None:
    src = Path(str(source or ""))
    if not src.is_file():
        raise FileNotFoundError(str(src))
    if src.resolve() == target.resolve():
        return
    if residue_name and target.suffix.lower() == ".mol2":
        target.write_text(
            _mcpb_mol2_with_residue_name(
                src.read_text(encoding="utf-8"),
                residue_name,
                atom_name_map,
            ),
            encoding="utf-8",
        )
        return
    shutil.copy2(src, target)


def _mcpb_mol2_with_residue_name(
    text: str,
    residue_name: str,
    atom_name_map: object,
) -> str:
    lines = text.splitlines()
    in_molecule = False
    in_atom = False
    molecule_name_done = False
    output: list[str] = []
    name_map = _string_map(atom_name_map)
    for line in lines:
        upper = line.upper()
        if upper.startswith("@<TRIPOS>MOLECULE"):
            in_molecule = True
            in_atom = False
            molecule_name_done = False
            output.append(line)
            continue
        if upper.startswith("@<TRIPOS>ATOM"):
            in_molecule = False
            in_atom = True
            output.append(line)
            continue
        if upper.startswith("@<TRIPOS>"):
            in_molecule = False
            in_atom = False
            output.append(line)
            continue
        if in_molecule and line.strip() and not molecule_name_done:
            output.append(residue_name)
            molecule_name_done = True
            continue
        if in_atom and line.strip():
            output.append(
                _mcpb_mol2_atom_line_with_residue_name(line, residue_name, name_map)
            )
            continue
        output.append(line)
    return "\n".join(output) + "\n"


def _mcpb_mol2_atom_line_with_residue_name(
    line: str,
    residue_name: str,
    atom_name_map: Mapping[str, str],
) -> str:
    fields = line.split()
    if len(fields) < 9:
        return line
    fields[1] = atom_name_map.get(fields[1], fields[1])
    fields[7] = residue_name
    return (
        f"{int(fields[0]):7d} {fields[1]:<8s} "
        f"{float(fields[2]):10.4f} {float(fields[3]):10.4f} {float(fields[4]):10.4f} "
        f"{fields[5]:<6s} {int(fields[6]):5d} {fields[7]:<8s} "
        f"{float(fields[8]):10.6f}"
    )


def _cofactor_mol2_atom_name_map(
    receptor_pdb: Path | None,
    residue: Mapping[str, str],
    sdf: Path,
    mol2: Path,
) -> dict[str, object]:
    pdb_names_by_sdf_order = _pdb_atom_names_by_sdf_order(
        receptor_pdb,
        residue,
        sdf,
        include_hydrogens=False,
    )
    mol2_names = _mol2_atom_names(mol2, include_hydrogens=False)
    all_atom_map = _cofactor_all_atom_name_map(receptor_pdb, residue, sdf, mol2)
    if not pdb_names_by_sdf_order or not mol2_names:
        return {
            "mol2_atom_name_map_status": "unavailable",
            "mol2_atom_name_map": {},
            **all_atom_map,
        }
    if len(pdb_names_by_sdf_order) != len(mol2_names):
        return {
            "mol2_atom_name_map_status": "atom_count_mismatch",
            "mol2_atom_name_map": {},
            **all_atom_map,
            "mol2_heavy_atom_count": len(mol2_names),
            "pdb_heavy_atom_count": len(pdb_names_by_sdf_order),
        }
    atom_name_map = dict(zip(mol2_names, pdb_names_by_sdf_order, strict=True))
    if isinstance(all_atom_map.get("mol2_all_atom_name_map"), Mapping):
        atom_name_map.update(_string_map(all_atom_map["mol2_all_atom_name_map"]))
    return {
        "mol2_atom_name_map_status": "mapped",
        "mol2_atom_name_map": atom_name_map,
        **all_atom_map,
        "mol2_heavy_atom_count": len(mol2_names),
        "pdb_heavy_atom_count": len(pdb_names_by_sdf_order),
    }


def _cofactor_all_atom_name_map(
    receptor_pdb: Path | None,
    residue: Mapping[str, str],
    sdf: Path,
    mol2: Path,
) -> dict[str, object]:
    pdb_names = _pdb_atom_names_by_sdf_order(
        receptor_pdb,
        residue,
        sdf,
        include_hydrogens=True,
    )
    mol2_names = _mol2_atom_names(mol2, include_hydrogens=True)
    if not pdb_names or not mol2_names:
        return {
            "mol2_all_atom_name_map_status": "unavailable",
            "mol2_all_atom_name_map": {},
        }
    if len(pdb_names) != len(mol2_names):
        return {
            "mol2_all_atom_name_map_status": "atom_count_mismatch",
            "mol2_all_atom_name_map": {},
            "mol2_all_atom_count": len(mol2_names),
            "pdb_all_atom_count": len(pdb_names),
        }
    return {
        "mol2_all_atom_name_map_status": "mapped",
        "mol2_all_atom_name_map": dict(zip(mol2_names, pdb_names, strict=True)),
        "mol2_all_atom_count": len(mol2_names),
        "pdb_all_atom_count": len(pdb_names),
    }


def _pdb_atom_names_by_sdf_order(
    receptor_pdb: Path | None,
    residue: Mapping[str, str],
    sdf: Path,
    *,
    include_hydrogens: bool,
) -> list[str]:
    if receptor_pdb is None or not receptor_pdb.exists():
        return []
    try:
        from rdkit import Chem
    except Exception:
        return []

    pdb_block = _pdb_residue_block(receptor_pdb, residue)
    if not pdb_block:
        return []
    pdb_mol = Chem.MolFromPDBBlock(
        pdb_block,
        sanitize=False,
        removeHs=False,
        proximityBonding=True,
    )
    sdf_mol = Chem.SDMolSupplier(str(sdf), sanitize=False, removeHs=False)[0]
    if pdb_mol is None or sdf_mol is None:
        return []
    return _pdb_names_for_matching_sdf_atoms(sdf_mol, pdb_mol, include_hydrogens)


def _pdb_residue_block(receptor_pdb: Path, residue: Mapping[str, str]) -> str:
    resname = str(residue.get("resname", "")).upper()
    chain = str(residue.get("chain", ""))
    resseq = str(residue.get("resseq", ""))
    lines: list[str] = []
    with receptor_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[17:20].strip().upper() != resname:
                continue
            if chain and line[21:22].strip() != chain:
                continue
            if resseq and line[22:26].strip() != resseq:
                continue
            lines.append(line)
    return "".join(lines) + "END\n" if lines else ""


def _pdb_names_for_matching_sdf_atoms(
    sdf_mol: Any,
    pdb_mol: Any,
    include_hydrogens: bool,
) -> list[str]:
    sdf_graph = _atom_graph(sdf_mol, include_hydrogens)
    pdb_graph = _atom_graph(pdb_mol, include_hydrogens)
    mapping = _match_heavy_atom_graphs(sdf_graph, pdb_graph)
    if not mapping:
        return []
    names = _pdb_atom_names(pdb_mol, include_hydrogens)
    return [names[mapping[idx]] for idx in range(len(sdf_graph["elements"]))]


def _atom_graph(mol: Any, include_hydrogens: bool) -> dict[str, Any]:
    old_ids = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if include_hydrogens or atom.GetSymbol().upper() not in {"H", "D"}
    ]
    old_to_new = {old: new for new, old in enumerate(old_ids)}
    neighbors: list[set[int]] = [set() for _ in old_ids]
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in old_to_new and end in old_to_new:
            neighbors[old_to_new[begin]].add(old_to_new[end])
            neighbors[old_to_new[end]].add(old_to_new[begin])
    return {
        "elements": [
            mol.GetAtomWithIdx(old).GetSymbol().upper() for old in old_ids
        ],
        "neighbors": neighbors,
    }


def _match_heavy_atom_graphs(
    source: Mapping[str, object],
    target: Mapping[str, object],
) -> dict[int, int]:
    source_elements = [str(item) for item in cast(Sequence[object], source["elements"])]
    target_elements = [str(item) for item in cast(Sequence[object], target["elements"])]
    source_neighbors = list(cast(Sequence[set[int]], source["neighbors"]))
    target_neighbors = list(cast(Sequence[set[int]], target["neighbors"]))
    if len(source_elements) != len(target_elements):
        return {}
    candidates = _graph_match_candidates(
        source_elements,
        target_elements,
        source_neighbors,
        target_neighbors,
    )
    if not all(candidates):
        return {}
    order = sorted(range(len(source_elements)), key=lambda idx: len(candidates[idx]))
    return _backtrack_graph_match(order, candidates, source_neighbors, target_neighbors)


def _graph_match_candidates(
    source_elements: Sequence[str],
    target_elements: Sequence[str],
    source_neighbors: Sequence[set[int]],
    target_neighbors: Sequence[set[int]],
) -> list[list[int]]:
    source_sigs = _graph_signatures(source_elements, source_neighbors)
    target_sigs = _graph_signatures(target_elements, target_neighbors)
    return [
        [
            target_idx
            for target_idx, target_sig in enumerate(target_sigs)
            if target_sig == source_sig
        ]
        for source_sig in source_sigs
    ]


def _graph_signatures(
    elements: Sequence[str],
    neighbors: Sequence[set[int]],
) -> list[tuple[str, int, tuple[str, ...]]]:
    return [
        (
            element,
            len(neighbors[idx]),
            tuple(sorted(elements[neighbor] for neighbor in neighbors[idx])),
        )
        for idx, element in enumerate(elements)
    ]


def _backtrack_graph_match(
    order: Sequence[int],
    candidates: Sequence[Sequence[int]],
    source_neighbors: Sequence[set[int]],
    target_neighbors: Sequence[set[int]],
) -> dict[int, int]:
    assigned: dict[int, int] = {}
    used: set[int] = set()

    def search(position: int) -> bool:
        if position == len(order):
            return True
        source_idx = order[position]
        for target_idx in candidates[source_idx]:
            if target_idx in used:
                continue
            if not _partial_graph_match_ok(
                source_idx,
                target_idx,
                assigned,
                source_neighbors,
                target_neighbors,
            ):
                continue
            assigned[source_idx] = target_idx
            used.add(target_idx)
            if search(position + 1):
                return True
            used.remove(target_idx)
            del assigned[source_idx]
        return False

    return dict(assigned) if search(0) else {}


def _partial_graph_match_ok(
    source_idx: int,
    target_idx: int,
    assigned: Mapping[int, int],
    source_neighbors: Sequence[set[int]],
    target_neighbors: Sequence[set[int]],
) -> bool:
    for assigned_source, assigned_target in assigned.items():
        source_bonded = assigned_source in source_neighbors[source_idx]
        target_bonded = assigned_target in target_neighbors[target_idx]
        if source_bonded != target_bonded:
            return False
    return True


def _pdb_atom_names(mol: Any, include_hydrogens: bool) -> list[str]:
    names: list[str] = []
    for atom in mol.GetAtoms():
        if not include_hydrogens and atom.GetSymbol().upper() in {"H", "D"}:
            continue
        info = atom.GetPDBResidueInfo()
        names.append(info.GetName().strip() if info else atom.GetSymbol())
    return names


def _mol2_atom_names(mol2: Path, *, include_hydrogens: bool) -> list[str]:
    names: list[str] = []
    for fields in _mol2_atom_fields(mol2):
        if not include_hydrogens and _mol2_atom_is_hydrogen(fields):
            continue
        names.append(fields[1])
    return names


def _mol2_atom_fields(mol2: Path) -> list[list[str]]:
    fields: list[list[str]] = []
    in_atom = False
    for line in mol2.read_text(encoding="utf-8").splitlines():
        upper = line.upper()
        if upper.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            continue
        if upper.startswith("@<TRIPOS>"):
            in_atom = False
        elif in_atom and line.strip():
            fields.append(line.split())
    return fields


def _mol2_atom_is_hydrogen(fields: Sequence[str]) -> bool:
    if len(fields) < 6:
        return False
    atom_name = fields[1].upper()
    atom_type = fields[5].upper()
    return atom_name.startswith(("H", "D")) or atom_type.startswith(("H", "D"))


def _string_map(raw: object) -> dict[str, str]:
    if not isinstance(raw, Mapping):
        return {}
    return {str(key): str(value) for key, value in raw.items() if str(key)}


def _site_nonprotein_resnames(site: Mapping[str, object]) -> set[str]:
    donors = site.get("donors", [])
    if not isinstance(donors, list):
        return set()
    resnames: set[str] = set()
    for donor in donors:
        if not isinstance(donor, Mapping):
            continue
        if str(donor.get("category", "")).lower() == "protein":
            continue
        resname = str(donor.get("resname", "")).upper()
        if resname:
            resnames.add(resname)
    return resnames


def _cofactor_files(
    entries: Sequence[Mapping[str, object]],
    key: str,
) -> list[str]:
    files: list[str] = []
    for entry in entries:
        path = str(entry.get(key, ""))
        if path:
            files.append(path)
    return files


def _parameterized_cofactor_resnames(cofactor_plan: Mapping[str, object]) -> list[str]:
    entries = cofactor_plan.get("entries", [])
    if not isinstance(entries, list):
        return []
    return sorted(
        {
            str(entry.get("resname", "")).upper()
            for entry in entries
            if isinstance(entry, Mapping) and bool(entry.get("ok"))
        }
    )


def _safe_site_id(site: Mapping[str, object], index: int) -> str:
    raw = str(site.get("id", "") or f"metal_{index}")
    safe = "".join(ch if ch.isalnum() else "_" for ch in raw).strip("_")
    return safe or f"metal_{index}"


def _donor_ids(donors: Sequence[object]) -> list[str]:
    ids: list[str] = []
    for donor in donors:
        if not isinstance(donor, Mapping):
            continue
        ids.append(
            ":".join(
                part
                for part in (
                    str(donor.get("resname", "")),
                    str(donor.get("chain", "")),
                    str(donor.get("resseq", "")),
                    str(donor.get("atom_name", "")),
                )
                if part
            )
        )
    return ids


def _mcpb_charge(site: Mapping[str, object]) -> str:
    charge = str(site.get("formal_charge_primary", "")).replace("+", "")
    return charge or "REVIEW"


def _mcpb_ion_info(site: Mapping[str, object], charge: str) -> str:
    resname = str(site.get("resname", "") or site.get("element", "")).upper()
    atom_name = str(site.get("atom_name", "") or site.get("element", "")).upper()
    element = str(site.get("element", "") or atom_name).upper()
    charge_token = charge if charge.lstrip("-").isdigit() else "REVIEW"
    return " ".join((resname[:3] or element[:2], atom_name[:4], element[:2], charge_token))


def _ion_id(site: Mapping[str, object], receptor_info: Mapping[str, object] | None = None) -> str:
    if receptor_info:
        serial_by_key = receptor_info.get("serial_by_key", {})
        if isinstance(serial_by_key, Mapping):
            mapped = str(serial_by_key.get(_site_atom_key(site), "")).strip()
            if mapped:
                return mapped
    for key in ("atom_serial", "serial", "metal_serial"):
        value = str(site.get(key, "")).strip()
        if value:
            return value
    return "REVIEW_METAL_ATOM_ID"


def _protonation_notes(site: Mapping[str, object]) -> list[str]:
    raw = site.get("metal_bound_protonation_recommendations", [])
    if not isinstance(raw, list):
        return []
    notes: list[str] = []
    for entry in raw:
        if not isinstance(entry, Mapping):
            continue
        residue = str(entry.get("residue", ""))
        recommendation = str(entry.get("recommendation", ""))
        suggested = str(entry.get("suggested_amber_resname", ""))
        note = f"{residue}:{recommendation}"
        if suggested:
            note += f":suggested={suggested}"
        notes.append(note)
    return notes


def _write_metal_pdb(
    work_dir: Path,
    site_id: str,
    site: Mapping[str, object],
) -> Path | None:
    coords = _coords(site)
    if coords is None:
        return None
    element = str(site.get("element", "") or "M").upper()[:2]
    ion_id = _ion_id(site)
    serial = _safe_int(ion_id) if ion_id.isdigit() else 1
    path = work_dir / f"{site_id}.pdb"
    x, y, z = coords
    line = (
        f"HETATM{serial:5d} {element:>4s} {element:>3s} A{serial % 9999:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
        "END\n"
    )
    path.write_text(line, encoding="utf-8")
    return path


def _write_mcpb_receptor(
    work_dir: Path,
    receptor_pdb: Path | None,
    preclean_plan: Mapping[str, object],
) -> dict[str, object]:
    if receptor_pdb is None or not receptor_pdb.exists():
        return {"path": "", "serial_by_key": {}, "preclean": preclean_plan}
    output = work_dir / "mcpb_receptor_precleaned.pdb"
    serial_by_key: dict[tuple[str, str, str, str], str] = {}
    overrides = _preclean_override_map(preclean_plan)
    drop_counts: Counter[str] = Counter()
    drop_residues: dict[str, Counter[str]] = {
        "drop_unparameterized_het": Counter(),
        "drop_unparameterized_cofactor": Counter(),
    }
    applied_overrides: dict[tuple[str, str, str], str] = {}
    residue_numbers: dict[tuple[str, str, str], int] = {}
    next_serial = 1
    next_resseq = 1
    lines: list[str] = []
    with receptor_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                lines.append(line)
                continue
            drop_reason = _mcpb_drop_reason(line, preclean_plan)
            if drop_reason:
                drop_counts[drop_reason] += 1
                if drop_reason in drop_residues:
                    drop_residues[drop_reason][line[17:20].strip().upper()] += 1
                continue
            serial = str(next_serial)
            atom_key = _line_atom_key(line)
            cleaned_line = _apply_residue_override(line, overrides)
            residue_identity = _line_residue_identity(line)
            if residue_identity not in residue_numbers:
                residue_numbers[residue_identity] = next_resseq
                next_resseq += 1
            cleaned_line = _renumber_residue(cleaned_line, residue_numbers[residue_identity])
            cleaned_key = _line_atom_key(cleaned_line)
            serial_by_key[atom_key] = serial
            serial_by_key[cleaned_key] = serial
            if atom_key[:3] in overrides:
                applied_overrides[atom_key[:3]] = overrides[atom_key[:3]]
            lines.append(f"{cleaned_line[:6]}{next_serial:5d}{cleaned_line[11:]}")
            next_serial += 1
    output.write_text("".join(lines), encoding="utf-8")
    oxt_repair = repair_terminal_oxt_geometry_in_pdb(
        output,
        audit_path=output.with_suffix(".oxt_repair_audit.json"),
    )
    preclean = dict(preclean_plan)
    preclean.update(
        {
            "status": "applied",
            "dropped_hydrogen_count": drop_counts["drop_hydrogen"],
            "dropped_terminal_oxt_count": 0,
            "terminal_oxt_policy": "preserve_and_repair_geometry",
            "terminal_oxt_repair": oxt_repair,
            "dropped_unparameterized_het_atom_count": drop_counts[
                "drop_unparameterized_het"
            ],
            "dropped_unparameterized_het_residue_counts": _format_counter(
                drop_residues["drop_unparameterized_het"]
            ),
            "dropped_unparameterized_cofactor_atom_count": drop_counts[
                "drop_unparameterized_cofactor"
            ],
            "dropped_unparameterized_cofactor_residue_counts": _format_counter(
                drop_residues["drop_unparameterized_cofactor"]
            ),
            "het_policy": (
                "keep standard residues, aliases.yaml metals/waters, and cofactors "
                "with generated MOL2/FRCMOD; drop remaining non-metal HET residues "
                "from heuristic MCPB staging"
            ),
            "renumbered_residue_count": len(residue_numbers),
            "residue_renumbering_policy": (
                "renumber MCPB staging-copy residues uniquely across chains because "
                "MCPB.py indexes model residues by residue number"
            ),
            "applied_residue_name_overrides": _serialize_residue_overrides(
                applied_overrides
            ),
            "applied_residue_name_override_count": len(applied_overrides),
        }
    )
    return {"path": str(output), "serial_by_key": serial_by_key, "preclean": preclean}


def _line_atom_key(line: str) -> tuple[str, str, str, str]:
    return (
        line[17:20].strip().upper(),
        line[21:22].strip(),
        line[22:26].strip(),
        line[12:16].strip().upper(),
    )


def _line_residue_identity(line: str) -> tuple[str, str, str]:
    return (line[21:22].strip(), line[22:26].strip(), line[26:27].strip())


def _site_atom_key(site: Mapping[str, object]) -> tuple[str, str, str, str]:
    return (
        str(site.get("resname", "")).upper(),
        str(site.get("chain", "")),
        str(site.get("resseq", "")),
        str(site.get("atom_name", "")).upper(),
    )


def _mcpb_preclean_plan(
    sites: Sequence[Mapping[str, object]],
    alias_sets: Mapping[str, set[str]],
    cofactor_plan: Mapping[str, object],
) -> dict[str, object]:
    suggestions: dict[tuple[str, str, str], set[str]] = {}
    for site in sites:
        raw = site.get("metal_bound_protonation_recommendations", [])
        if not isinstance(raw, list):
            continue
        for entry in raw:
            if not isinstance(entry, Mapping):
                continue
            residue_key = _recommendation_residue_key(entry)
            suggested = str(entry.get("suggested_amber_resname", "")).upper()
            if not residue_key or not suggested:
                continue
            suggestions.setdefault(residue_key, set()).add(suggested)

    overrides: dict[tuple[str, str, str], str] = {}
    conflicts: list[dict[str, object]] = []
    for residue_key, names in sorted(suggestions.items()):
        if len(names) == 1:
            overrides[residue_key] = next(iter(names))
        else:
            conflicts.append(
                {
                    "residue": _format_residue_key(residue_key),
                    "suggested_amber_resnames": sorted(names),
                }
            )

    return {
        "status": "planned",
        "alias_metal_count": len(alias_sets.get("metals", set())),
        "alias_water_count": len(alias_sets.get("waters", set())),
        "alias_cofactor_count": len(alias_sets.get("cofactors", set())),
        "hydrogen_policy": "drop_all_hydrogens_from_mcpb_staging_copy",
        "residue_name_policy": (
            "apply unambiguous local metal-bound Amber residue-name heuristics "
            "only in the MCPB staging copy"
        ),
        "planned_residue_name_overrides": _serialize_residue_overrides(overrides),
        "planned_residue_name_override_count": len(overrides),
        "residue_name_conflicts": conflicts,
        "residue_name_conflict_count": len(conflicts),
        "parameterized_cofactor_resnames": _parameterized_cofactor_resnames(
            cofactor_plan
        ),
        "canonical_cofactor_resnames": sorted(alias_sets.get("cofactors", set())),
        "canonical_metal_resnames": sorted(alias_sets.get("metals", set())),
        "canonical_water_resnames": sorted(alias_sets.get("waters", set())),
        "protonation_scope": "local_metal_bound_heuristic_not_global_hbond_optimization",
    }


def _recommendation_residue_key(
    recommendation: Mapping[str, object],
) -> tuple[str, str, str] | None:
    residue = str(recommendation.get("residue", ""))
    parts = residue.split(":")
    if len(parts) < 3:
        return None
    return (parts[0].upper(), parts[1], parts[2])


def _preclean_override_map(
    preclean_plan: Mapping[str, object],
) -> dict[tuple[str, str, str], str]:
    raw = preclean_plan.get("planned_residue_name_overrides", [])
    if not isinstance(raw, list):
        return {}
    overrides: dict[tuple[str, str, str], str] = {}
    for entry in raw:
        if not isinstance(entry, Mapping):
            continue
        resname = str(entry.get("resname", "")).upper()
        chain = str(entry.get("chain", ""))
        resseq = str(entry.get("resseq", ""))
        suggested = str(entry.get("suggested_amber_resname", "")).upper()
        if resname and resseq and suggested:
            overrides[(resname, chain, resseq)] = suggested
    return overrides


def _apply_residue_override(
    line: str,
    overrides: Mapping[tuple[str, str, str], str],
) -> str:
    residue_key = (line[17:20].strip().upper(), line[21:22].strip(), line[22:26].strip())
    replacement = overrides.get(residue_key, "")
    if not replacement:
        return line
    return f"{line[:17]}{replacement[:3]:>3s}{line[20:]}"


def _renumber_residue(line: str, resseq: int) -> str:
    return f"{line[:22]}{resseq:4d} {line[27:]}"


def _is_hydrogen_record(line: str) -> bool:
    element = line[76:78].strip().upper()
    atom_name = line[12:16].strip().upper()
    if element in {"H", "D"}:
        return True
    return atom_name.startswith(("H", "D"))


def _mcpb_drop_reason(line: str, preclean_plan: Mapping[str, object]) -> str:
    if _is_hydrogen_record(line):
        return "drop_hydrogen"
    action = _mcpb_het_action(line, preclean_plan)
    return "" if action == "keep" else action


def _mcpb_het_action(line: str, preclean_plan: Mapping[str, object]) -> str:
    if not line.startswith("HETATM"):
        return "keep"
    resname = line[17:20].strip().upper()
    element = line[76:78].strip().upper()
    waters = _plan_token_set(preclean_plan, "canonical_water_resnames", WATER_RESNAMES)
    metals = _plan_token_set(preclean_plan, "canonical_metal_resnames", METAL_RESNAMES)
    cofactors = _plan_token_set(preclean_plan, "canonical_cofactor_resnames", set())
    parameterized = _plan_token_set(
        preclean_plan, "parameterized_cofactor_resnames", set()
    )
    if resname in STANDARD_RESNAMES or resname in waters:
        return "keep"
    if resname in metals or element in metals:
        return "keep"
    if resname in parameterized:
        return "keep"
    if resname in cofactors:
        return "drop_unparameterized_cofactor"
    return "drop_unparameterized_het"


def _plan_token_set(
    plan: Mapping[str, object],
    key: str,
    fallback: set[str],
) -> set[str]:
    raw = plan.get(key, [])
    if not isinstance(raw, list):
        return set(fallback)
    return {str(token).strip().upper() for token in raw if str(token).strip()}


def _serialize_residue_overrides(
    overrides: Mapping[tuple[str, str, str], str],
) -> list[dict[str, str]]:
    serialized: list[dict[str, str]] = []
    for (resname, chain, resseq), suggested in sorted(overrides.items()):
        serialized.append(
            {
                "residue": _format_residue_key((resname, chain, resseq)),
                "resname": resname,
                "chain": chain,
                "resseq": resseq,
                "suggested_amber_resname": suggested,
            }
        )
    return serialized


def _format_residue_key(residue_key: tuple[str, str, str]) -> str:
    resname, chain, resseq = residue_key
    return ":".join(part for part in (resname, chain, resseq) if part)


def _preclean_notes(
    receptor_info: Mapping[str, object],
    site: Mapping[str, object],
) -> list[str]:
    preclean = receptor_info.get("preclean", {})
    if not isinstance(preclean, Mapping):
        return []
    notes = [
        str(preclean.get("hydrogen_policy", "")),
        str(preclean.get("het_policy", "")),
        str(preclean.get("residue_renumbering_policy", "")),
        str(preclean.get("protonation_scope", "")),
    ]
    notes.extend(_preclean_override_notes(preclean, site))
    notes.extend(_preclean_count_notes(preclean))
    return [note for note in notes if note]


def _preclean_override_notes(
    preclean: Mapping[str, object],
    site: Mapping[str, object],
) -> list[str]:
    site_residues = _site_recommendation_residues(site)
    raw_overrides = preclean.get("applied_residue_name_overrides", [])
    if not isinstance(raw_overrides, list):
        return []
    return [
        f"{entry.get('residue', '')}->Amber {entry.get('suggested_amber_resname', '')}"
        for entry in raw_overrides
        if isinstance(entry, Mapping)
        and str(entry.get("residue", "")) in site_residues
        and str(entry.get("suggested_amber_resname", ""))
    ]


def _preclean_count_notes(preclean: Mapping[str, object]) -> list[str]:
    count_fields = (
        ("dropped_hydrogens", "dropped_hydrogen_count"),
        ("dropped_terminal_oxt", "dropped_terminal_oxt_count"),
        (
            "dropped_unparameterized_cofactor_atoms",
            "dropped_unparameterized_cofactor_atom_count",
        ),
        (
            "dropped_unparameterized_het_atoms",
            "dropped_unparameterized_het_atom_count",
        ),
    )
    return [
        f"{label}={count}"
        for label, field in count_fields
        if (count := _safe_int(preclean.get(field, 0)))
    ]


def _site_recommendation_residues(site: Mapping[str, object]) -> set[str]:
    residues: set[str] = set()
    raw = site.get("metal_bound_protonation_recommendations", [])
    if not isinstance(raw, list):
        return residues
    for entry in raw:
        if not isinstance(entry, Mapping):
            continue
        residue_key = _recommendation_residue_key(entry)
        if residue_key:
            residues.add(_format_residue_key(residue_key))
    return residues


def _summarize_preclean(preclean: object) -> dict[str, object]:
    if not isinstance(preclean, Mapping):
        return {
            "mcpb_preclean_status": "",
            "mcpb_dropped_hydrogen_count": 0,
            "mcpb_dropped_terminal_oxt_count": 0,
            "mcpb_dropped_unparameterized_cofactor_atom_count": 0,
            "mcpb_dropped_unparameterized_het_atom_count": 0,
            "mcpb_renumbered_residue_count": 0,
            "mcpb_residue_name_override_count": 0,
            "mcpb_residue_name_conflict_count": 0,
        }
    return {
        "mcpb_preclean_status": str(preclean.get("status", "")),
        "mcpb_dropped_hydrogen_count": _safe_int(
            preclean.get("dropped_hydrogen_count", 0)
        ),
        "mcpb_dropped_terminal_oxt_count": _safe_int(
            preclean.get("dropped_terminal_oxt_count", 0)
        ),
        "mcpb_dropped_unparameterized_cofactor_atom_count": _safe_int(
            preclean.get("dropped_unparameterized_cofactor_atom_count", 0)
        ),
        "mcpb_dropped_unparameterized_het_atom_count": _safe_int(
            preclean.get("dropped_unparameterized_het_atom_count", 0)
        ),
        "mcpb_renumbered_residue_count": _safe_int(
            preclean.get("renumbered_residue_count", 0)
        ),
        "mcpb_residue_name_override_count": _safe_int(
            preclean.get("applied_residue_name_override_count", 0)
        ),
        "mcpb_residue_name_conflict_count": _safe_int(
            preclean.get("residue_name_conflict_count", 0)
        ),
    }


def _summarize_cofactor_plan(plan: object) -> dict[str, object]:
    if not isinstance(plan, Mapping):
        return {
            "mcpb_cofactor_parameterization_status": "",
            "mcpb_cofactor_residue_count": 0,
            "mcpb_parameterized_cofactor_count": 0,
            "mcpb_unresolved_cofactor_count": 0,
        }
    return {
        "mcpb_cofactor_parameterization_status": str(plan.get("status", "")),
        "mcpb_cofactor_residue_count": _safe_int(
            plan.get("cofactor_residue_count", 0)
        ),
        "mcpb_parameterized_cofactor_count": _safe_int(
            plan.get("parameterized_count", 0)
        ),
        "mcpb_unresolved_cofactor_count": _safe_int(plan.get("unresolved_count", 0)),
        "mcpb_cofactor_atom_name_map_status_counts": _cofactor_atom_map_counts(plan),
    }


def _cofactor_atom_map_counts(plan: Mapping[str, object]) -> str:
    entries = plan.get("entries", [])
    if not isinstance(entries, list):
        return ""
    return _format_counter(
        Counter(
            str(entry.get("mol2_atom_name_map_status", "unknown"))
            for entry in entries
            if isinstance(entry, Mapping)
        )
    )


def _coords(site: Mapping[str, object]) -> tuple[float, float, float] | None:
    raw = site.get("coords") or site.get("metal_coords")
    if not isinstance(raw, (list, tuple)) or len(raw) != 3:
        return None
    try:
        x, y, z = raw
        return float(x), float(y), float(z)
    except Exception:
        return None


def _format_counter(counter: Counter[str]) -> str:
    return ",".join(f"{key}:{counter[key]}" for key in sorted(counter) if key)


def _safe_int(value: object) -> int:
    try:
        return int(float(str(value or 0)))
    except Exception:
        return 0
