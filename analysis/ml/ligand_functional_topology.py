"""Structure-only RDKit functional, branch, and topological-shape descriptors."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Callable


DescriptorValue = float | int
DescriptorCalculator = Callable[[Any], DescriptorValue]


@dataclass(frozen=True)
class MolecularDescriptor:
    column: str
    definition: str
    calculate: DescriptorCalculator


@dataclass(frozen=True)
class MolecularDescriptorGroup:
    name: str
    descriptors: tuple[MolecularDescriptor, ...]

    @property
    def columns(self) -> tuple[str, ...]:
        return tuple(descriptor.column for descriptor in self.descriptors)

    def calculate(self, mol: Any) -> dict[str, DescriptorValue]:
        return {
            descriptor.column: descriptor.calculate(mol)
            for descriptor in self.descriptors
        }


def _fragment_count(fragment_name: str) -> DescriptorCalculator:
    def calculate(mol: Any) -> int:
        from rdkit.Chem import Fragments

        fragment = getattr(Fragments, fragment_name)
        return int(fragment(mol, countUnique=True))

    return calculate


@lru_cache(maxsize=None)
def _smarts_query(smarts: str) -> Any:
    from rdkit import Chem

    query = Chem.MolFromSmarts(smarts)
    if query is None:
        raise ValueError(f"invalid descriptor SMARTS: {smarts}")
    return query


def _smarts_count(smarts: str) -> DescriptorCalculator:
    def calculate(mol: Any) -> int:
        matches = mol.GetSubstructMatches(_smarts_query(smarts), uniquify=True)
        return int(len(matches))

    return calculate


def _graph_descriptor(descriptor_name: str) -> DescriptorCalculator:
    def calculate(mol: Any) -> float:
        from rdkit.Chem import GraphDescriptors

        descriptor = getattr(GraphDescriptors, descriptor_name)
        return float(descriptor(mol))

    return calculate


# These partition RDKit's fr_COO pattern by the oxygen's explicit graph state.
_CARBOXYLIC_ACID = _smarts_count("[#6][CX3](=O)[O;H1;+0]")
_CARBOXYLATE = _smarts_count("[#6][CX3](=O)[O;-1]")

# Derived from RDKit's functional-group hierarchy. Classification uses carbon
# substitution and single bonds so explicit protonation does not change the type.
_PRIMARY_AMINE = _smarts_count("[N;D1;$(N-!@[#6]);!$(N-C=[O,N,S])]")
_SECONDARY_AMINE = _smarts_count("[N;D2;$(N(-[#6])-[#6]);!$(N-C=[O,N,S])]")
_TERTIARY_AMINE = _smarts_count("[N;D3;$(N(-[#6])(-[#6])-[#6]);!$(N-C=[O,N,S])]")


def _amine_count(mol: Any) -> int:
    return int(_PRIMARY_AMINE(mol) + _SECONDARY_AMINE(mol) + _TERTIARY_AMINE(mol))


def _heavy_degree(atom: Any) -> int:
    return sum(neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors())


def _heavy_atom_branch_points(mol: Any) -> int:
    return sum(
        atom.GetAtomicNum() > 1 and _heavy_degree(atom) >= 3 for atom in mol.GetAtoms()
    )


def _terminal_heavy_atoms(mol: Any) -> int:
    return sum(
        atom.GetAtomicNum() > 1 and _heavy_degree(atom) == 1 for atom in mol.GetAtoms()
    )


def _ring_nonring_attachment_bonds(mol: Any) -> int:
    count = 0
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if begin.GetAtomicNum() <= 1 or end.GetAtomicNum() <= 1:
            continue
        count += begin.IsInRing() != end.IsInRing()
    return count


def _nonring_heavy_atoms(mol: Any) -> int:
    return sum(
        atom.GetAtomicNum() > 1 and not atom.IsInRing() for atom in mol.GetAtoms()
    )


def _outside_bemis_murcko_scaffold(mol: Any) -> int:
    from rdkit.Chem.Scaffolds import MurckoScaffold

    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return int(mol.GetNumHeavyAtoms() - scaffold.GetNumHeavyAtoms())


def _bemis_murcko_side_chain_fragments(mol: Any) -> int:
    """Count disconnected substituent fragments outside a non-empty scaffold."""
    from rdkit.Chem.Scaffolds import MurckoScaffold

    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms() == 0:
        return 0
    scaffold_atoms = set(mol.GetSubstructMatch(scaffold))
    if len(scaffold_atoms) != scaffold.GetNumAtoms():
        return 0
    outside_atoms = {
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() not in scaffold_atoms
    }
    components = 0
    while outside_atoms:
        components += 1
        frontier = [outside_atoms.pop()]
        while frontier:
            atom = mol.GetAtomWithIdx(frontier.pop())
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in outside_atoms:
                    outside_atoms.remove(neighbor_idx)
                    frontier.append(neighbor_idx)
    return components


FUNCTIONAL_COUNT_GROUP = MolecularDescriptorGroup(
    name="functional_counts",
    descriptors=(
        MolecularDescriptor(
            "rdkit_aldehyde_count",
            "RDKit Fragments.fr_aldehyde unique-match count",
            _fragment_count("fr_aldehyde"),
        ),
        MolecularDescriptor(
            "rdkit_ester_count",
            "RDKit Fragments.fr_ester unique-match count",
            _fragment_count("fr_ester"),
        ),
        MolecularDescriptor(
            "rdkit_ketone_count",
            "RDKit Fragments.fr_ketone unique-match count",
            _fragment_count("fr_ketone"),
        ),
        MolecularDescriptor(
            "rdkit_amide_count",
            "RDKit Fragments.fr_amide unique-match count",
            _fragment_count("fr_amide"),
        ),
        MolecularDescriptor(
            "rdkit_thiol_count",
            "RDKit Fragments.fr_SH unique-match count",
            _fragment_count("fr_SH"),
        ),
        MolecularDescriptor(
            "rdkit_carboxylic_acid_count",
            "neutral [#6][CX3](=O)[O;H1;+0] partition of RDKit fr_COO",
            _CARBOXYLIC_ACID,
        ),
        MolecularDescriptor(
            "rdkit_carboxylate_count",
            "anionic [#6][CX3](=O)[O;-1] partition of RDKit fr_COO",
            _CARBOXYLATE,
        ),
        MolecularDescriptor(
            "rdkit_amine_count",
            "sum of primary, secondary, and tertiary carbon-substituted amines",
            _amine_count,
        ),
        MolecularDescriptor(
            "rdkit_primary_amine_count",
            "amine N with one single-bonded carbon substituent",
            _PRIMARY_AMINE,
        ),
        MolecularDescriptor(
            "rdkit_secondary_amine_count",
            "amine N with two single-bonded carbon substituents",
            _SECONDARY_AMINE,
        ),
        MolecularDescriptor(
            "rdkit_tertiary_amine_count",
            "amine N with three single-bonded carbon substituents",
            _TERTIARY_AMINE,
        ),
        MolecularDescriptor(
            "rdkit_benzene_ring_count",
            "RDKit Fragments.fr_benzene c1ccccc1 unique-match count",
            _fragment_count("fr_benzene"),
        ),
    ),
)

BRANCH_COUNT_GROUP = MolecularDescriptorGroup(
    name="branch_counts",
    descriptors=(
        MolecularDescriptor(
            "rdkit_heavy_atom_branch_point_count",
            "heavy atoms with at least three heavy-atom neighbors",
            _heavy_atom_branch_points,
        ),
        MolecularDescriptor(
            "rdkit_terminal_heavy_atom_count",
            "heavy atoms with exactly one heavy-atom neighbor",
            _terminal_heavy_atoms,
        ),
        MolecularDescriptor(
            "rdkit_ring_nonring_attachment_bond_count",
            "heavy-atom bonds with exactly one ring-atom endpoint",
            _ring_nonring_attachment_bonds,
        ),
        MolecularDescriptor(
            "rdkit_nonring_heavy_atom_count",
            "heavy atoms for which RDKit Atom.IsInRing is false",
            _nonring_heavy_atoms,
        ),
        MolecularDescriptor(
            "rdkit_bemis_murcko_outside_heavy_atom_count",
            "molecule heavy atoms minus Bemis-Murcko scaffold heavy atoms",
            _outside_bemis_murcko_scaffold,
        ),
        MolecularDescriptor(
            "rdkit_bemis_murcko_side_chain_fragment_count",
            "connected substituent fragments outside a non-empty Bemis-Murcko scaffold",
            _bemis_murcko_side_chain_fragments,
        ),
    ),
)

TOPOLOGICAL_SHAPE_GROUP = MolecularDescriptorGroup(
    name="topological_shape",
    descriptors=(
        MolecularDescriptor(
            "rdkit_bertz_ct",
            "RDKit GraphDescriptors.BertzCT molecular complexity index",
            _graph_descriptor("BertzCT"),
        ),
        MolecularDescriptor(
            "rdkit_hall_kier_alpha",
            "RDKit GraphDescriptors.HallKierAlpha correction factor",
            _graph_descriptor("HallKierAlpha"),
        ),
        MolecularDescriptor(
            "rdkit_kappa1",
            "RDKit GraphDescriptors.Kappa1 shape index",
            _graph_descriptor("Kappa1"),
        ),
        MolecularDescriptor(
            "rdkit_kappa2",
            "RDKit GraphDescriptors.Kappa2 shape index",
            _graph_descriptor("Kappa2"),
        ),
        MolecularDescriptor(
            "rdkit_kappa3",
            "RDKit GraphDescriptors.Kappa3 shape index",
            _graph_descriptor("Kappa3"),
        ),
    ),
)

LIGAND_STRUCTURE_DESCRIPTOR_GROUPS = (
    FUNCTIONAL_COUNT_GROUP,
    BRANCH_COUNT_GROUP,
    TOPOLOGICAL_SHAPE_GROUP,
)
FUNCTIONAL_COUNT_COLUMNS = FUNCTIONAL_COUNT_GROUP.columns
BRANCH_COUNT_COLUMNS = BRANCH_COUNT_GROUP.columns
TOPOLOGICAL_SHAPE_COLUMNS = TOPOLOGICAL_SHAPE_GROUP.columns
FUNCTIONAL_TOPOLOGY_COLUMNS = tuple(
    column for group in LIGAND_STRUCTURE_DESCRIPTOR_GROUPS for column in group.columns
)
FUNCTIONAL_TOPOLOGY_DEFINITIONS = {
    descriptor.column: descriptor.definition
    for group in LIGAND_STRUCTURE_DESCRIPTOR_GROUPS
    for descriptor in group.descriptors
}


def functional_count_record(mol: Any) -> dict[str, DescriptorValue]:
    """Return registered functional-group counts for a molecule."""
    return FUNCTIONAL_COUNT_GROUP.calculate(mol)


def branch_count_record(mol: Any) -> dict[str, DescriptorValue]:
    """Return exact heavy-atom branch and scaffold counts for a molecule."""
    return BRANCH_COUNT_GROUP.calculate(mol)


def topological_shape_record(mol: Any) -> dict[str, DescriptorValue]:
    """Return established RDKit graph-complexity and shape descriptors."""
    return TOPOLOGICAL_SHAPE_GROUP.calculate(mol)


def functional_topology_record(mol: Any) -> dict[str, DescriptorValue]:
    """Return every registered functional, branch, and shape descriptor."""
    record: dict[str, DescriptorValue] = {}
    for group in LIGAND_STRUCTURE_DESCRIPTOR_GROUPS:
        record.update(group.calculate(mol))
    return record
