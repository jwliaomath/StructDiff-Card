from __future__ import annotations

from collections import defaultdict
from itertools import product
from math import acos, degrees, floor, sqrt

import numpy as np

from .models import (
    ContactChange,
    ContactSummary,
    MotionDecomposition,
    ResidueMapping,
    SpatialPatch,
)


def _coordinates(row: ResidueMapping, *, mobile: bool) -> np.ndarray | None:
    values = (
        (row.mobile_aligned_x, row.mobile_aligned_y, row.mobile_aligned_z)
        if mobile
        else (row.reference_x, row.reference_y, row.reference_z)
    )
    if any(value is None for value in values):
        return None
    return np.asarray(values, dtype=float)


def _nearby_pairs(coordinates: dict[int, np.ndarray], cutoff: float) -> set[tuple[int, int]]:
    """Return pairs within cutoff using a small spatial hash instead of an N² matrix."""

    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    pairs: set[tuple[int, int]] = set()
    offsets = list(product((-1, 0, 1), repeat=3))
    for index in sorted(coordinates):
        point = coordinates[index]
        cell = tuple(floor(float(value) / cutoff) for value in point)
        for offset in offsets:
            neighbor = tuple(cell[axis] + offset[axis] for axis in range(3))
            for other in cells.get(neighbor, []):
                if float(np.linalg.norm(point - coordinates[other])) <= cutoff:
                    pairs.add((other, index))
        cells[cell].append(index)
    return pairs


def summarize_contact_changes(
    rows: list[ResidueMapping],
    *,
    contact_cutoff: float = 8.0,
    transition_margin: float = 1.0,
    sequence_exclusion: int = 2,
    max_changes: int = 500,
) -> tuple[ContactSummary, list[ContactChange]]:
    """Find stable representative-atom contacts gained or lost relative to reference.

    A contact is at most ``contact_cutoff`` in one structure and at least the cutoff plus
    ``transition_margin`` in the other. The margin suppresses changes caused by tiny
    fluctuations at the threshold. Immediate same-chain sequence neighbors are excluded.
    """

    mapped = [
        row
        for row in rows
        if row.status == "matched"
        and _coordinates(row, mobile=True) is not None
        and _coordinates(row, mobile=False) is not None
    ]
    mobile_coordinates = {
        index: coordinate
        for index, row in enumerate(mapped)
        if (coordinate := _coordinates(row, mobile=True)) is not None
    }
    reference_coordinates = {
        index: coordinate
        for index, row in enumerate(mapped)
        if (coordinate := _coordinates(row, mobile=False)) is not None
    }
    candidate_pairs = _nearby_pairs(mobile_coordinates, contact_cutoff) | _nearby_pairs(
        reference_coordinates, contact_cutoff
    )

    local_positions: dict[int, int] = {}
    chain_cursors: dict[tuple[str | None, str | None], int] = defaultdict(int)
    for index, row in enumerate(mapped):
        chain_pair = (row.mobile_chain, row.reference_chain)
        local_positions[index] = chain_cursors[chain_pair]
        chain_cursors[chain_pair] += 1

    retained = 0
    gained = 0
    lost = 0
    inter_chain_gained = 0
    inter_chain_lost = 0
    raw_changes: list[tuple[float, bool, str, int, int, float, float]] = []
    far_threshold = contact_cutoff + transition_margin
    for first_index, second_index in candidate_pairs:
        first = mapped[first_index]
        second = mapped[second_index]
        same_chain_pair = (
            first.mobile_chain == second.mobile_chain
            and first.reference_chain == second.reference_chain
        )
        if same_chain_pair and abs(local_positions[first_index] - local_positions[second_index]) <= (
            sequence_exclusion
        ):
            continue

        mobile_distance = float(
            np.linalg.norm(mobile_coordinates[first_index] - mobile_coordinates[second_index])
        )
        reference_distance = float(
            np.linalg.norm(reference_coordinates[first_index] - reference_coordinates[second_index])
        )
        mobile_contact = mobile_distance <= contact_cutoff
        reference_contact = reference_distance <= contact_cutoff
        inter_chain = (
            first.mobile_chain != second.mobile_chain
            or first.reference_chain != second.reference_chain
        )
        if mobile_contact and reference_contact:
            retained += 1
            continue
        if mobile_contact and reference_distance >= far_threshold:
            change = "gained"
            gained += 1
            inter_chain_gained += int(inter_chain)
        elif reference_contact and mobile_distance >= far_threshold:
            change = "lost"
            lost += 1
            inter_chain_lost += int(inter_chain)
        else:
            continue
        raw_changes.append(
            (
                abs(mobile_distance - reference_distance),
                inter_chain,
                change,
                first_index,
                second_index,
                mobile_distance,
                reference_distance,
            )
        )

    raw_changes.sort(key=lambda item: (item[1], item[0]), reverse=True)
    changes: list[ContactChange] = []
    for rank, item in enumerate(raw_changes[:max_changes], start=1):
        _, inter_chain, change, first_index, second_index, mobile_distance, reference_distance = item
        first = mapped[first_index]
        second = mapped[second_index]
        changes.append(
            ContactChange(
                rank=rank,
                change=change,
                alignment_index_1=first.alignment_index,
                alignment_index_2=second.alignment_index,
                mobile_chain_1=first.mobile_chain or "?",
                mobile_residue_1=first.mobile_residue or "?",
                mobile_chain_2=second.mobile_chain or "?",
                mobile_residue_2=second.mobile_residue or "?",
                reference_chain_1=first.reference_chain or "?",
                reference_residue_1=first.reference_residue or "?",
                reference_chain_2=second.reference_chain or "?",
                reference_residue_2=second.reference_residue or "?",
                mobile_distance=mobile_distance,
                reference_distance=reference_distance,
                distance_change=mobile_distance - reference_distance,
                inter_chain=inter_chain,
            )
        )

    summary = ContactSummary(
        mapped_residues=len(mapped),
        contact_cutoff=contact_cutoff,
        transition_margin=transition_margin,
        retained_contacts=retained,
        gained_contacts=gained,
        lost_contacts=lost,
        inter_chain_gained_contacts=inter_chain_gained,
        inter_chain_lost_contacts=inter_chain_lost,
    )
    return summary, changes


def decompose_chain_motion(rows: list[ResidueMapping]) -> list[MotionDecomposition]:
    """Separate each chain's global residual into local rigid and internal components."""

    grouped: dict[tuple[str, str], list[ResidueMapping]] = defaultdict(list)
    for row in rows:
        if (
            row.status == "matched"
            and row.mobile_chain is not None
            and row.reference_chain is not None
            and _coordinates(row, mobile=True) is not None
            and _coordinates(row, mobile=False) is not None
        ):
            grouped[(row.mobile_chain, row.reference_chain)].append(row)

    results: list[MotionDecomposition] = []
    for (mobile_chain, reference_chain), chain_rows in grouped.items():
        mobile = np.vstack([_coordinates(row, mobile=True) for row in chain_rows])
        reference = np.vstack([_coordinates(row, mobile=False) for row in chain_rows])
        global_rmsd = float(sqrt(np.mean(np.sum((mobile - reference) ** 2, axis=1))))
        mobile_centroid = np.mean(mobile, axis=0)
        reference_centroid = np.mean(reference, axis=0)
        centered_mobile = mobile - mobile_centroid
        centered_reference = reference - reference_centroid
        covariance = centered_mobile.T @ centered_reference
        left, _, right_transpose = np.linalg.svd(covariance)
        rotation = right_transpose.T @ left.T
        if np.linalg.det(rotation) < 0:
            right_transpose[-1, :] *= -1
            rotation = right_transpose.T @ left.T
        fitted = centered_mobile @ rotation.T + reference_centroid
        internal_rmsd = float(sqrt(np.mean(np.sum((fitted - reference) ** 2, axis=1))))
        rigid_body_rmsd = sqrt(max(global_rmsd**2 - internal_rmsd**2, 0.0))
        rigid_fraction = rigid_body_rmsd**2 / global_rmsd**2 if global_rmsd else 0.0
        centroid_shift = float(np.linalg.norm(mobile_centroid - reference_centroid))
        angle_argument = float(np.clip((np.trace(rotation) - 1.0) / 2.0, -1.0, 1.0))
        rotation_degrees = degrees(acos(angle_argument))

        if len(chain_rows) < 3:
            classification = "insufficient anchors"
        elif global_rmsd < 1.0:
            classification = "well aligned"
        elif rigid_fraction >= 0.67:
            classification = "mostly rigid-body"
        elif rigid_fraction <= 0.33:
            classification = "mostly internal"
        else:
            classification = "mixed"
        results.append(
            MotionDecomposition(
                mobile_chain=mobile_chain,
                reference_chain=reference_chain,
                matched_residues=len(chain_rows),
                global_rmsd=global_rmsd,
                internal_rmsd=internal_rmsd,
                rigid_body_rmsd=rigid_body_rmsd,
                rigid_body_fraction=rigid_fraction,
                centroid_shift=centroid_shift,
                correction_rotation_degrees=rotation_degrees,
                classification=classification,
            )
        )
    return results


def cluster_spatial_patches(
    rows: list[ResidueMapping],
    *,
    threshold: float,
    neighbor_cutoff: float = 10.0,
    max_patches: int = 20,
) -> list[SpatialPatch]:
    """Group high-residual residues by 3D proximity in the reference frame."""

    hot_rows = [
        row
        for row in rows
        if row.displacement is not None
        and row.displacement >= threshold
        and _coordinates(row, mobile=False) is not None
    ]
    coordinates = {
        index: coordinate
        for index, row in enumerate(hot_rows)
        if (coordinate := _coordinates(row, mobile=False)) is not None
    }
    adjacency: dict[int, set[int]] = {index: set() for index in coordinates}
    for first, second in _nearby_pairs(coordinates, neighbor_cutoff):
        adjacency[first].add(second)
        adjacency[second].add(first)

    components: list[list[int]] = []
    unseen = set(coordinates)
    while unseen:
        seed = unseen.pop()
        component = [seed]
        stack = [seed]
        while stack:
            current = stack.pop()
            neighbors = adjacency[current] & unseen
            unseen.difference_update(neighbors)
            stack.extend(neighbors)
            component.extend(neighbors)
        components.append(component)

    components.sort(
        key=lambda component: (
            max(float(hot_rows[index].displacement or 0.0) for index in component),
            len(component),
        ),
        reverse=True,
    )
    patches: list[SpatialPatch] = []
    for rank, component in enumerate(components[:max_patches], start=1):
        component_rows = [hot_rows[index] for index in component]
        points = np.vstack([coordinates[index] for index in component])
        centroid = np.mean(points, axis=0)
        values = [float(row.displacement or 0.0) for row in component_rows]
        mobile_chains = sorted({row.mobile_chain or "?" for row in component_rows})
        reference_chains = sorted({row.reference_chain or "?" for row in component_rows})
        patches.append(
            SpatialPatch(
                rank=rank,
                member_count=len(component_rows),
                mobile_chains=", ".join(mobile_chains),
                reference_chains=", ".join(reference_chains),
                mobile_residues=", ".join(
                    f"{row.mobile_chain}:{row.mobile_residue}" for row in component_rows
                ),
                reference_residues=", ".join(
                    f"{row.reference_chain}:{row.reference_residue}" for row in component_rows
                ),
                centroid_x=float(centroid[0]),
                centroid_y=float(centroid[1]),
                centroid_z=float(centroid[2]),
                mean_displacement=sum(values) / len(values),
                peak_displacement=max(values),
                spatial_radius=max(
                    (float(np.linalg.norm(point - centroid)) for point in points), default=0.0
                ),
            )
        )
    return patches
