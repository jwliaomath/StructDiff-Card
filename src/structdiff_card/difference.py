from __future__ import annotations

from collections import defaultdict
from statistics import median

from .models import DifferenceRegion, DifferenceSummary, ResidueMapping


def _percentile(values: list[float], quantile: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    position = (len(ordered) - 1) * quantile
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def summarize_differences(
    rows: list[ResidueMapping],
    *,
    absolute_threshold: float = 2.0,
    context_residues: int = 2,
    max_regions: int = 8,
) -> tuple[DifferenceSummary, list[DifferenceRegion]]:
    """Summarize coherent residual regions after the one global alignment transform.

    Hot residues use the larger of 2 Å and the 90th percentile. Adjacent hot
    residues separated by at most one cooler residue are grouped, then expanded
    by two aligned residues on each side to provide structural context.
    """

    matched = [row for row in rows if row.displacement is not None]
    values = [float(row.displacement) for row in matched if row.displacement is not None]
    p90 = _percentile(values, 0.90)
    threshold = max(absolute_threshold, p90)
    hot_count = sum(value >= threshold for value in values)
    summary = DifferenceSummary(
        matched_residues=len(values),
        median_displacement=median(values) if values else 0.0,
        p90_displacement=p90,
        max_displacement=max(values, default=0.0),
        hotspot_threshold=threshold,
        residues_above_threshold=hot_count,
        fraction_above_threshold=hot_count / len(values) if values else 0.0,
    )

    grouped: dict[tuple[str, str], list[ResidueMapping]] = defaultdict(list)
    for row in matched:
        if row.mobile_chain is not None and row.reference_chain is not None:
            grouped[(row.mobile_chain, row.reference_chain)].append(row)

    candidates: list[tuple[float, float, tuple[str, str], list[ResidueMapping], int]] = []
    for chain_pair, chain_rows in grouped.items():
        hot_indices = [
            index
            for index, row in enumerate(chain_rows)
            if row.displacement is not None and row.displacement >= threshold
        ]
        if not hot_indices:
            continue

        segments: list[list[int]] = [[hot_indices[0]]]
        for index in hot_indices[1:]:
            previous = segments[-1][-1]
            matched_gap = index - previous
            alignment_gap = (
                chain_rows[index].alignment_index - chain_rows[previous].alignment_index
            )
            if matched_gap <= 2 and alignment_gap == matched_gap:
                segments[-1].append(index)
            else:
                segments.append([index])

        expanded: list[tuple[int, int, int]] = []
        for segment in segments:
            start = segment[0]
            end = segment[-1]
            for _ in range(context_residues):
                if (
                    start > 0
                    and chain_rows[start].alignment_index
                    - chain_rows[start - 1].alignment_index
                    == 1
                ):
                    start -= 1
                if (
                    end < len(chain_rows) - 1
                    and chain_rows[end + 1].alignment_index
                    - chain_rows[end].alignment_index
                    == 1
                ):
                    end += 1
            high_count = len(segment)
            merge_with_previous = False
            if expanded:
                previous_end = expanded[-1][1]
                merge_with_previous = start <= previous_end or (
                    start == previous_end + 1
                    and chain_rows[start].alignment_index
                    - chain_rows[previous_end].alignment_index
                    == 1
                )
            if merge_with_previous:
                old_start, old_end, old_count = expanded[-1]
                expanded[-1] = (old_start, max(old_end, end), old_count + high_count)
            else:
                expanded.append((start, end, high_count))

        for start, end, high_count in expanded:
            region_rows = chain_rows[start : end + 1]
            region_values = [float(row.displacement) for row in region_rows if row.displacement]
            candidates.append(
                (
                    max(region_values),
                    sum(region_values) / len(region_values),
                    chain_pair,
                    region_rows,
                    high_count,
                )
            )

    candidates.sort(key=lambda item: (item[0], item[1]), reverse=True)
    regions: list[DifferenceRegion] = []
    for rank, (_, mean_value, chain_pair, region_rows, high_count) in enumerate(
        candidates[:max_regions], start=1
    ):
        peak = max(region_rows, key=lambda row: row.displacement or 0.0)
        first, last = region_rows[0], region_rows[-1]
        regions.append(
            DifferenceRegion(
                rank=rank,
                mobile_chain=chain_pair[0],
                reference_chain=chain_pair[1],
                mobile_start=first.mobile_residue or "?",
                mobile_end=last.mobile_residue or "?",
                reference_start=first.reference_residue or "?",
                reference_end=last.reference_residue or "?",
                peak_mobile_residue=peak.mobile_residue or "?",
                peak_reference_residue=peak.reference_residue or "?",
                mean_displacement=mean_value,
                peak_displacement=float(peak.displacement or 0.0),
                residue_count=len(region_rows),
                residues_above_threshold=high_count,
            )
        )

    return summary, regions
