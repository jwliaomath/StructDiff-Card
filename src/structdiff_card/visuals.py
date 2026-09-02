from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import plotly.graph_objects as go

from .models import AlignmentResult
from .structure import read_structure

_MOBILE_COLORS = ["#31d6c6", "#23a7b7", "#6be6d1", "#178d9b", "#91f2df"]
_REFERENCE_COLORS = ["#ff8a5c", "#e86f4a", "#ffb07f", "#d65f3d", "#ffc39f"]
_DIFFERENCE_COLORSCALE = [
    [0.00, "#2166ac"],
    [0.20, "#67a9cf"],
    [0.40, "#d1e5f0"],
    [0.60, "#fddbc7"],
    [0.80, "#ef8a62"],
    [1.00, "#b2182b"],
]


def overlay_figure(
    mobile_path: str | Path,
    reference_path: str | Path,
    result: AlignmentResult,
    *,
    all_models: bool = False,
) -> go.Figure:
    mobile = read_structure(mobile_path, all_models=all_models)
    reference = read_structure(reference_path, all_models=all_models)
    figure = go.Figure()
    displacement_lookup = {
        (row.mobile_chain, row.mobile_residue): row.displacement
        for row in result.residue_mapping
        if row.mobile_chain is not None and row.mobile_residue is not None
    }
    show_scale = True

    for chain, residues in reference.chains.items():
        xyz = [residue.xyz for residue in residues]
        figure.add_trace(
            go.Scatter3d(
                x=[item[0] for item in xyz],
                y=[item[1] for item in xyz],
                z=[item[2] for item in xyz],
                mode="lines",
                name=f"Reference · {chain}",
                line={"color": "#d8e2e6", "width": 7},
                opacity=0.30,
                hovertemplate=f"Reference chain {chain}<extra></extra>",
            )
        )
    for chain, residues in mobile.chains.items():
        xyz = [result.transform.apply(residue.xyz) for residue in residues]
        if chain in result.unmatched_mobile_chains:
            figure.add_trace(
                go.Scatter3d(
                    x=[item[0] for item in xyz],
                    y=[item[1] for item in xyz],
                    z=[item[2] for item in xyz],
                    mode="lines",
                    name=f"Unmatched mobile · {chain}",
                    line={"color": "#b77ce8", "width": 5},
                    opacity=0.42,
                    hovertemplate=f"Mobile chain {chain} · unmatched<extra></extra>",
                )
            )
            continue

        values = [displacement_lookup.get((chain, residue.residue_id)) for residue in residues]
        numeric_values = [value if value is not None else 0.0 for value in values]
        labels = [
            f"{chain}:{residue.residue_id} · "
            + (f"{value:.2f} Å" if value is not None else "not mapped")
            for residue, value in zip(residues, values, strict=True)
        ]
        figure.add_trace(
            go.Scatter3d(
                x=[item[0] for item in xyz],
                y=[item[1] for item in xyz],
                z=[item[2] for item in xyz],
                mode="lines+markers",
                name=f"Difference map · {chain}",
                line={"color": "rgba(213,226,229,0.48)", "width": 4},
                marker={
                    "size": 2.2,
                    "color": numeric_values,
                    "colorscale": _DIFFERENCE_COLORSCALE,
                    "cmin": 0,
                    "cmax": 5,
                    "showscale": show_scale,
                    "colorbar": {
                        "title": {
                            "text": "Residual<br>distance (Å)",
                            "font": {"color": "#c7d7dd"},
                        },
                        "tickvals": [0, 1, 2, 3, 4, 5],
                        "ticktext": ["0", "1", "2", "3", "4", "≥5"],
                        "thickness": 13,
                        "len": 0.55,
                        "x": 0.98,
                        "tickfont": {"color": "#c7d7dd"},
                    },
                },
                customdata=labels,
                hovertemplate="%{customdata}<extra></extra>",
            )
        )
        show_scale = False

    mobile_lookup = {
        (chain, residue.residue_id): residue
        for chain, residues in mobile.chains.items()
        for residue in residues
    }
    reference_lookup = {
        (chain, residue.residue_id): residue
        for chain, residues in reference.chains.items()
        for residue in residues
    }
    threshold = result.difference_summary.hotspot_threshold if result.difference_summary else 2.0
    vector_rows = sorted(
        (
            row
            for row in result.residue_mapping
            if row.displacement is not None and row.displacement >= threshold
        ),
        key=lambda row: row.displacement or 0.0,
        reverse=True,
    )[:20]
    vector_x: list[float | None] = []
    vector_y: list[float | None] = []
    vector_z: list[float | None] = []
    for row in vector_rows:
        mobile_residue = mobile_lookup.get((row.mobile_chain, row.mobile_residue))
        reference_residue = reference_lookup.get((row.reference_chain, row.reference_residue))
        if mobile_residue is None or reference_residue is None:
            continue
        moved = result.transform.apply(mobile_residue.xyz)
        vector_x.extend([reference_residue.xyz[0], moved[0], None])
        vector_y.extend([reference_residue.xyz[1], moved[1], None])
        vector_z.extend([reference_residue.xyz[2], moved[2], None])
    if vector_x:
        figure.add_trace(
            go.Scatter3d(
                x=vector_x,
                y=vector_y,
                z=vector_z,
                mode="lines",
                name=f"Difference vectors · ≥{threshold:.2f} Å",
                line={"color": "#ffcf56", "width": 5},
                hoverinfo="skip",
            )
        )

    visible_patches = [patch for patch in result.spatial_patches if patch.member_count >= 2][:8]
    if visible_patches:
        figure.add_trace(
            go.Scatter3d(
                x=[patch.centroid_x for patch in visible_patches],
                y=[patch.centroid_y for patch in visible_patches],
                z=[patch.centroid_z for patch in visible_patches],
                mode="markers+text",
                name="3D difference patches",
                marker={
                    "size": 5,
                    "symbol": "diamond",
                    "color": "#ffcf56",
                    "line": {"color": "#07131b", "width": 1},
                },
                text=[f"P{patch.rank}" for patch in visible_patches],
                textposition="top center",
                textfont={"color": "#ffcf56", "size": 11},
                customdata=[
                    f"Patch P{patch.rank} · {patch.member_count} residues · "
                    f"peak {patch.peak_displacement:.2f} Å"
                    for patch in visible_patches
                ],
                hovertemplate="%{customdata}<extra></extra>",
            )
        )

    figure.update_layout(
        margin={"l": 0, "r": 0, "t": 16, "b": 0},
        paper_bgcolor="#07131b",
        plot_bgcolor="#07131b",
        legend={"font": {"color": "#c7d7dd"}, "bgcolor": "rgba(0,0,0,0)"},
        scene={
            "bgcolor": "#07131b",
            "aspectmode": "data",
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "zaxis": {"visible": False},
            "camera": {"eye": {"x": 1.45, "y": 1.45, "z": 0.9}},
        },
        height=620,
    )
    return figure


def displacement_figure(result: AlignmentResult) -> go.Figure:
    grouped: dict[str, list] = defaultdict(list)
    for row in result.residue_mapping:
        if row.displacement is not None:
            grouped[f"{row.mobile_chain}→{row.reference_chain}"].append(row)

    figure = go.Figure()
    cursor = 0
    show_scale = True
    position_lookup: dict[tuple[str | None, str | None], tuple[int, float]] = {}
    for index, (chain_pair, rows) in enumerate(grouped.items()):
        x = list(range(cursor, cursor + len(rows)))
        labels = [
            f"{row.mobile_chain}:{row.mobile_residue} ↔ "
            f"{row.reference_chain}:{row.reference_residue}"
            for row in rows
        ]
        values = [float(row.displacement or 0.0) for row in rows]
        for point_x, row in zip(x, rows, strict=True):
            position_lookup[(row.mobile_chain, row.mobile_residue)] = (
                point_x,
                float(row.displacement or 0.0),
            )
        figure.add_trace(
            go.Scatter(
                x=x,
                y=values,
                mode="lines+markers",
                name=f"{chain_pair} · residues",
                line={"width": 1.3, "color": "#bdcbc9"},
                marker={
                    "size": 5,
                    "color": values,
                    "colorscale": _DIFFERENCE_COLORSCALE,
                    "cmin": 0,
                    "cmax": 5,
                    "showscale": show_scale,
                    "colorbar": {"title": "Residual (Å)", "thickness": 12},
                },
                customdata=labels,
                hovertemplate="%{customdata}<br>%{y:.2f} Å<extra></extra>",
            )
        )
        show_scale = False
        half_window = 3
        rolling = []
        for row_index in range(len(values)):
            start = max(0, row_index - half_window)
            end = min(len(values), row_index + half_window + 1)
            rolling.append(sum(values[start:end]) / (end - start))
        figure.add_trace(
            go.Scatter(
                x=x,
                y=rolling,
                mode="lines",
                name=f"{chain_pair} · 7-residue mean",
                line={"width": 3, "color": _MOBILE_COLORS[index % len(_MOBILE_COLORS)]},
                hovertemplate="7-residue mean: %{y:.2f} Å<extra></extra>",
            )
        )
        cursor += len(rows) + 3

    threshold = result.difference_summary.hotspot_threshold if result.difference_summary else 2.0
    figure.add_hline(
        y=threshold,
        line_dash="dot",
        line_color="#ff8a5c",
        annotation_text=f"hotspot threshold · {threshold:.2f} Å",
        annotation_font_color="#a94728",
    )
    hotspot_x = []
    hotspot_y = []
    hotspot_labels = []
    hotspot_ranks = []
    for region in result.difference_regions:
        point = position_lookup.get((region.mobile_chain, region.peak_mobile_residue))
        if point is None:
            continue
        hotspot_x.append(point[0])
        hotspot_y.append(point[1])
        hotspot_ranks.append(str(region.rank))
        hotspot_labels.append(
            f"Region {region.rank}: {region.mobile_chain}:{region.mobile_start}–"
            f"{region.mobile_end} · peak {region.peak_displacement:.2f} Å"
        )
    if hotspot_x:
        figure.add_trace(
            go.Scatter(
                x=hotspot_x,
                y=hotspot_y,
                mode="markers+text",
                name="Difference regions",
                marker={"size": 11, "symbol": "diamond", "color": "#b2182b"},
                text=hotspot_ranks,
                textposition="top center",
                customdata=hotspot_labels,
                hovertemplate="%{customdata}<extra></extra>",
            )
        )
    figure.update_layout(
        height=360,
        margin={"l": 52, "r": 20, "t": 20, "b": 48},
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        xaxis_title="Matched residue pairs (chain mapping order)",
        yaxis_title="Displacement (Å)",
        hovermode="closest",
        legend={"orientation": "h", "y": 1.04},
    )
    figure.update_xaxes(showgrid=False)
    figure.update_yaxes(gridcolor="#e1e9e7", zeroline=False)
    return figure


def contact_change_figure(result: AlignmentResult) -> go.Figure:
    figure = go.Figure()
    styles = {
        "gained": {"color": "#0d9488", "symbol": "circle", "label": "Gained in mobile"},
        "lost": {"color": "#e76f51", "symbol": "x", "label": "Lost from mobile"},
    }
    for change, style in styles.items():
        rows = [row for row in result.contact_changes if row.change == change]
        if not rows:
            continue
        labels = [
            f"{style['label']}<br>mobile {row.mobile_chain_1}:{row.mobile_residue_1} — "
            f"{row.mobile_chain_2}:{row.mobile_residue_2}: {row.mobile_distance:.2f} Å<br>"
            f"reference {row.reference_chain_1}:{row.reference_residue_1} — "
            f"{row.reference_chain_2}:{row.reference_residue_2}: "
            f"{row.reference_distance:.2f} Å"
            for row in rows
        ]
        figure.add_trace(
            go.Scatter(
                x=[row.alignment_index_1 for row in rows],
                y=[row.alignment_index_2 for row in rows],
                mode="markers",
                name=str(style["label"]),
                marker={
                    "size": [10 if row.inter_chain else 7 for row in rows],
                    "symbol": style["symbol"],
                    "color": style["color"],
                    "line": {"width": 1, "color": "#ffffff"},
                },
                customdata=labels,
                hovertemplate="%{customdata}<extra></extra>",
            )
        )
    if not result.contact_changes:
        figure.add_annotation(
            text="No stable gained/lost contacts at the selected 8 Å definition",
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
            font={"color": "#61747c"},
        )
    figure.update_layout(
        height=430,
        margin={"l": 62, "r": 24, "t": 28, "b": 58},
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        xaxis_title="Mapped residue index 1",
        yaxis_title="Mapped residue index 2",
        legend={"orientation": "h", "y": 1.08},
    )
    figure.update_xaxes(gridcolor="#edf2f1", zeroline=False)
    figure.update_yaxes(gridcolor="#edf2f1", zeroline=False, scaleanchor="x", scaleratio=1)
    return figure


def motion_decomposition_figure(result: AlignmentResult) -> go.Figure:
    labels = [
        f"{row.mobile_chain}→{row.reference_chain}" for row in result.motion_decomposition
    ]
    figure = go.Figure()
    figure.add_trace(
        go.Bar(
            x=labels,
            y=[row.rigid_body_rmsd for row in result.motion_decomposition],
            name="Rigid-body component",
            marker_color="#287fb8",
            customdata=[row.classification for row in result.motion_decomposition],
            hovertemplate="%{x}<br>Rigid-body: %{y:.2f} Å<br>%{customdata}<extra></extra>",
        )
    )
    figure.add_trace(
        go.Bar(
            x=labels,
            y=[row.internal_rmsd for row in result.motion_decomposition],
            name="Internal deformation",
            marker_color="#e76f51",
            customdata=[row.classification for row in result.motion_decomposition],
            hovertemplate="%{x}<br>Internal: %{y:.2f} Å<br>%{customdata}<extra></extra>",
        )
    )
    figure.update_layout(
        barmode="group",
        height=360,
        margin={"l": 52, "r": 20, "t": 26, "b": 46},
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        xaxis_title="Chain mapping",
        yaxis_title="RMSD component (Å)",
        legend={"orientation": "h", "y": 1.08},
    )
    figure.update_xaxes(showgrid=False)
    figure.update_yaxes(gridcolor="#e1e9e7", zeroline=False)
    return figure


def gap_track_html(result: AlignmentResult) -> str:
    sections: list[str] = []
    for alignment in result.chain_alignments:
        cells: list[str] = []
        for mobile, reference in zip(
            alignment.mobile_sequence, alignment.reference_sequence, strict=True
        ):
            if mobile == "-":
                css = "only-reference"
                label = reference
                title = "Only in reference"
            elif reference == "-":
                css = "only-mobile"
                label = mobile
                title = "Only in mobile"
            else:
                css = "matched"
                label = mobile
                title = "Matched"
            cells.append(f'<span class="track-cell {css}" title="{title}">{label}</span>')
        sections.append(
            f'<div class="track-row"><div class="track-label">'
            f'{alignment.mobile_chain} → {alignment.reference_chain}</div>'
            f'<div class="track-sequence">{"".join(cells)}</div></div>'
        )

    unmatched = []
    for chain in result.unmatched_mobile_chains:
        unmatched.append(f'<span class="unmatched-chip">Mobile {chain} · unmatched</span>')
    for chain in result.unmatched_reference_chains:
        unmatched.append(f'<span class="unmatched-chip">Reference {chain} · unmatched</span>')
    return (
        '<div class="track-wrap">'
        + "".join(sections)
        + f'<div class="unmatched-row">{"".join(unmatched) or "No unmatched chains"}</div>'
        + "</div>"
    )
