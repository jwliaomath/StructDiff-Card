from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class Transform:
    """Rigid transform applied to Structure 1 (mobile) to overlay Structure 2."""

    rotation: tuple[tuple[float, float, float], ...]
    translation: tuple[float, float, float]

    def apply(self, xyz: tuple[float, float, float]) -> tuple[float, float, float]:
        x, y, z = xyz
        r = self.rotation
        t = self.translation
        return (
            t[0] + r[0][0] * x + r[0][1] * y + r[0][2] * z,
            t[1] + r[1][0] * x + r[1][1] * y + r[1][2] * z,
            t[2] + r[2][0] * x + r[2][1] * y + r[2][2] * z,
        )


@dataclass(frozen=True)
class ScoreSummary:
    mobile_length: int
    reference_length: int
    aligned_length: int
    rmsd: float
    sequence_identity: float
    tm_mobile: float
    tm_reference: float


@dataclass(frozen=True)
class ChainAlignment:
    mobile_chain: str
    reference_chain: str
    mobile_length: int
    reference_length: int
    aligned_length: int
    rmsd: float
    sequence_identity: float
    tm_mobile: float
    tm_reference: float
    mobile_sequence: str
    marker: str
    reference_sequence: str
    coverage_mobile: float = 0.0
    coverage_reference: float = 0.0


@dataclass(frozen=True)
class ResidueMapping:
    alignment_index: int
    status: str
    mobile_chain: str | None
    mobile_residue: str | None
    mobile_resname: str | None
    mobile_one_letter: str | None
    reference_chain: str | None
    reference_residue: str | None
    reference_resname: str | None
    reference_one_letter: str | None
    displacement: float | None
    delta_x: float | None = None
    delta_y: float | None = None
    delta_z: float | None = None
    mobile_aligned_x: float | None = None
    mobile_aligned_y: float | None = None
    mobile_aligned_z: float | None = None
    reference_x: float | None = None
    reference_y: float | None = None
    reference_z: float | None = None


@dataclass(frozen=True)
class DifferenceSummary:
    matched_residues: int
    median_displacement: float
    p90_displacement: float
    max_displacement: float
    hotspot_threshold: float
    residues_above_threshold: int
    fraction_above_threshold: float


@dataclass(frozen=True)
class DifferenceRegion:
    rank: int
    mobile_chain: str
    reference_chain: str
    mobile_start: str
    mobile_end: str
    reference_start: str
    reference_end: str
    peak_mobile_residue: str
    peak_reference_residue: str
    mean_displacement: float
    peak_displacement: float
    residue_count: int
    residues_above_threshold: int


@dataclass(frozen=True)
class ContactSummary:
    mapped_residues: int
    contact_cutoff: float
    transition_margin: float
    retained_contacts: int
    gained_contacts: int
    lost_contacts: int
    inter_chain_gained_contacts: int
    inter_chain_lost_contacts: int


@dataclass(frozen=True)
class ContactChange:
    rank: int
    change: str
    alignment_index_1: int
    alignment_index_2: int
    mobile_chain_1: str
    mobile_residue_1: str
    mobile_chain_2: str
    mobile_residue_2: str
    reference_chain_1: str
    reference_residue_1: str
    reference_chain_2: str
    reference_residue_2: str
    mobile_distance: float
    reference_distance: float
    distance_change: float
    inter_chain: bool


@dataclass(frozen=True)
class MotionDecomposition:
    mobile_chain: str
    reference_chain: str
    matched_residues: int
    global_rmsd: float
    internal_rmsd: float
    rigid_body_rmsd: float
    rigid_body_fraction: float
    centroid_shift: float
    correction_rotation_degrees: float
    classification: str


@dataclass(frozen=True)
class SpatialPatch:
    rank: int
    member_count: int
    mobile_chains: str
    reference_chains: str
    mobile_residues: str
    reference_residues: str
    centroid_x: float
    centroid_y: float
    centroid_z: float
    mean_displacement: float
    peak_displacement: float
    spatial_radius: float


@dataclass
class AlignmentResult:
    project_version: str
    engine: str
    engine_version: str
    mode: str
    assembly_selection: str
    command: list[str]
    mobile_name: str
    reference_name: str
    mobile_sha256: str
    reference_sha256: str
    summary: ScoreSummary
    transform: Transform
    chain_alignments: list[ChainAlignment]
    residue_mapping: list[ResidueMapping]
    difference_summary: DifferenceSummary | None = None
    difference_regions: list[DifferenceRegion] = field(default_factory=list)
    contact_summary: ContactSummary | None = None
    contact_changes: list[ContactChange] = field(default_factory=list)
    motion_decomposition: list[MotionDecomposition] = field(default_factory=list)
    spatial_patches: list[SpatialPatch] = field(default_factory=list)
    unmatched_mobile_chains: list[str] = field(default_factory=list)
    unmatched_reference_chains: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    raw_stdout: str = ""

    def to_dict(self, *, include_raw: bool = False) -> dict[str, Any]:
        payload = asdict(self)
        if not include_raw:
            payload.pop("raw_stdout", None)
        return payload
