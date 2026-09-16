from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from structdiff_card.analysis import compare_structures
from structdiff_card.export import write_bundle
from structdiff_card.settings import DEFAULT_TIMEOUT_SECONDS, configured_timeout_seconds
from structdiff_card.usalign import USAlignTimeoutError
from structdiff_card.visuals import (
    contact_change_figure,
    displacement_figure,
    gap_track_html,
    motion_decomposition_figure,
    overlay_figure,
)

ROOT = Path(__file__).resolve().parent
EXAMPLES = ROOT / "examples" / "input"

timeout_config_warning: str | None = None
try:
    configured_timeout = configured_timeout_seconds()
except ValueError as error:
    configured_timeout = DEFAULT_TIMEOUT_SECONDS
    timeout_config_warning = f"{error}; using {DEFAULT_TIMEOUT_SECONDS} seconds instead."

st.set_page_config(
    page_title="StructDiff Card",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
<style>
  :root { --ink: #07131b; --teal: #0d9488; --coral: #e76f51; --mist: #edf5f3; }
  .stApp { background: #f7faf9; color: #10242d; }
  [data-testid="stSidebar"] { background: #07131b; color: #edf7f5; }
  [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2,
  [data-testid="stSidebar"] h3, [data-testid="stSidebar"] label,
  [data-testid="stSidebar"] [data-testid="stCaptionContainer"],
  [data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p { color: #edf7f5 !important; }
  [data-testid="stSidebar"] [data-testid="stFileUploader"] section,
  [data-testid="stSidebar"] [data-testid="stFileUploaderFile"] {
    background: #eef5f3; border-color: #c9d9d5;
  }
  [data-testid="stSidebar"] [data-testid="stFileUploader"] section *,
  [data-testid="stSidebar"] [data-testid="stFileUploaderFile"] * { color: #17303a !important; }
  [data-testid="stSidebar"] [data-testid="stFileUploaderFile"] svg,
  [data-testid="stSidebar"] [data-testid="stFileUploader"] section svg { fill: #17303a; color: #17303a; }
  [data-testid="stSidebar"] .stTextInput input {
    background: #f7faf9; color: #10242d; border-color: #c9d9d5;
  }
  [data-testid="stSidebar"] .stTextInput input::placeholder { color: #63777f; opacity: 1; }
  [data-testid="stSidebar"] .stButton button[kind="secondary"] {
    background: #f7faf9; color: #17303a; border-color: #c9d9d5;
  }
  [data-testid="stSidebar"] .stButton button[kind="secondary"] p { color: #17303a !important; }
  [data-testid="stSidebar"] .stButton button[kind="primary"] p { color: #ffffff !important; }
  .hero-kicker { letter-spacing: .16em; text-transform: uppercase; color: #0d9488; font-weight: 700; font-size: .76rem; }
  .hero-title { font-family: Georgia, serif; font-size: clamp(2.3rem, 5vw, 4.8rem); line-height: .96; margin: .3rem 0 .8rem; }
  .hero-copy { color: #52666f; max-width: 780px; font-size: 1.05rem; }
  .metric { border-top: 3px solid #0d9488; padding-top: .8rem; }
  .metric-label { color: #61747c; text-transform: uppercase; letter-spacing: .1em; font-size: .72rem; }
  .metric-value { font-family: Georgia, serif; font-size: 2.2rem; color: #10242d; }
  .track-wrap { background: white; border: 1px solid #d8e4e1; padding: 1rem; overflow-x: auto; }
  .track-row { margin-bottom: 1rem; min-width: max-content; }
  .track-label { font: 700 .8rem/1.4 Arial; color: #415963; margin-bottom: .35rem; }
  .track-sequence { display: flex; gap: 1px; }
  .track-cell { width: 14px; height: 22px; display: inline-grid; place-items: center; font: 10px/1 monospace; }
  .track-cell.matched { background: #d5f3ed; color: #0d655e; }
  .track-cell.only-mobile { background: #ffe0d2; color: #a94728; }
  .track-cell.only-reference { background: #dce9ff; color: #31568d; }
  .unmatched-chip { display: inline-block; border-left: 4px solid #e76f51; background: #fff0e9; padding: .45rem .7rem; margin-right: .5rem; color: #7d321e; font-size: .8rem; }
  .provenance { font-family: monospace; font-size: .78rem; background: #eaf1ef; padding: .7rem; overflow-wrap: anywhere; }
  .publication-card svg { width: 100%; height: auto; display: block; }
  .difference-scale { height: 10px; background: linear-gradient(90deg, #2166ac, #67a9cf 20%, #d1e5f0 40%, #fddbc7 60%, #ef8a62 80%, #b2182b); margin: .35rem 0; }
  .difference-scale-labels { display: flex; justify-content: space-between; color: #61747c; font-size: .72rem; }
</style>
""",
    unsafe_allow_html=True,
)


def parse_manual_mapping(text: str) -> list[tuple[str, str]] | None:
    if not text.strip():
        return None
    pairs = []
    for item in text.split(","):
        if ":" not in item:
            raise ValueError("Manual mapping must look like A:C,B:D")
        mobile, reference = item.split(":", 1)
        pairs.append((mobile.strip(), reference.strip()))
    return pairs


def run_analysis(
    mobile_bytes: bytes,
    mobile_name: str,
    reference_bytes: bytes,
    reference_name: str,
    assembly: str,
    mapping_text: str,
    timeout_seconds: int | None = 1800,
) -> None:
    with tempfile.TemporaryDirectory(prefix="structdiff-ui-") as temp_dir:
        temp = Path(temp_dir)
        # Uploads commonly share a basename (e.g. two training runs' 10.pdb).
        # Separate input roles so writing the reference cannot replace mobile.
        mobile_dir = temp / "mobile_input"
        reference_dir = temp / "reference_input"
        mobile_dir.mkdir()
        reference_dir.mkdir()
        mobile_path = mobile_dir / Path(mobile_name).name
        reference_path = reference_dir / Path(reference_name).name
        mobile_path.write_bytes(mobile_bytes)
        reference_path.write_bytes(reference_bytes)
        result = compare_structures(
            mobile_path,
            reference_path,
            assembly_selection=assembly,
            manual_chain_mapping=parse_manual_mapping(mapping_text),
            timeout=timeout_seconds,
        )
        output = temp / "result"
        files = write_bundle(result, mobile_path, reference_path, output)
        st.session_state["analysis"] = {
            "result": result,
            "overlay": overlay_figure(
                mobile_path,
                reference_path,
                result,
                all_models=assembly == "biological",
            ),
            "displacement": displacement_figure(result),
            "contacts": contact_change_figure(result),
            "motion": motion_decomposition_figure(result),
            "files": {name: path.read_bytes() for name, path in files.items()},
        }


with st.sidebar:
    st.markdown("## Run locally")
    st.caption("Files stay in this process and are deleted after analysis.")
    mobile_upload = st.file_uploader("Structure 1 · mobile", type=["pdb", "cif", "mmcif"])
    reference_upload = st.file_uploader("Structure 2 · reference", type=["pdb", "cif", "mmcif"])
    assembly = st.radio(
        "Coordinate scope",
        options=["asymmetric", "biological"],
        format_func=lambda value: "First model / asymmetric unit"
        if value == "asymmetric"
        else "All models / biological assembly",
    )
    mapping_text = st.text_input(
        "Optional chain mapping",
        placeholder="C:A,D:B",
        help="Leave blank for US-align automatic chain assignment.",
    )
    with st.expander("Advanced settings"):
        if timeout_config_warning:
            st.warning(timeout_config_warning)
        timeout_seconds = st.number_input(
            "US-align timeout (seconds)",
            min_value=0,
            value=configured_timeout,
            step=300,
            help=(
                "Default is 1800 seconds (30 minutes). Set to 0 to allow US-align to run "
                "without a time limit."
            ),
        )
    analyze = st.button("Analyze structures", type="primary", width="stretch")
    example = st.button("Load 2HHB → 1HHO example", width="stretch")

    if analyze:
        if not mobile_upload or not reference_upload:
            st.error("Upload both structures first.")
        else:
            try:
                with st.spinner("Running US-align and building residue map…"):
                    run_analysis(
                        mobile_upload.getvalue(),
                        mobile_upload.name,
                        reference_upload.getvalue(),
                        reference_upload.name,
                        assembly,
                        mapping_text,
                        int(timeout_seconds),
                    )
            except USAlignTimeoutError as error:
                st.error(str(error))
            except Exception as error:  # noqa: BLE001 - UI boundary must display engine errors
                st.exception(error)
    if example:
        try:
            with st.spinner("Running the bundled assembly-mismatch example…"):
                run_analysis(
                    (EXAMPLES / "2HHB.pdb").read_bytes(),
                    "2HHB.pdb",
                    (EXAMPLES / "1HHO.pdb").read_bytes(),
                    "1HHO.pdb",
                    "asymmetric",
                    "",
                    int(timeout_seconds),
                )
        except USAlignTimeoutError as error:
            st.error(str(error))
        except Exception as error:  # noqa: BLE001 - UI boundary must display engine errors
            st.exception(error)

st.markdown('<div class="hero-kicker">Local-first structural comparison</div>', unsafe_allow_html=True)
st.markdown('<div class="hero-title">See where the assembly moves.</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="hero-copy">StructDiff Card turns a US-align superposition into an auditable '
    'chain map, contact-change map, motion decomposition, spatial hotspots, and publication-ready export. '
    'Unequal residue and chain counts remain visible instead of being hidden by one RMSD.</div>',
    unsafe_allow_html=True,
)

analysis = st.session_state.get("analysis")
if not analysis:
    st.info("Use the bundled example for an immediate result, or upload two PDB/mmCIF files.")
    st.stop()

result = analysis["result"]
summary = result.summary
metric_columns = st.columns(4)
metrics = [
    ("TM-score · reference", f"{summary.tm_reference:.3f}"),
    ("Global RMSD", f"{summary.rmsd:.2f} Å"),
    ("Aligned residues", f"{summary.aligned_length} / {summary.reference_length}"),
    ("Mapped chains", f"{len(result.chain_alignments)}"),
]
for column, (label, value) in zip(metric_columns, metrics, strict=True):
    column.markdown(
        f'<div class="metric"><div class="metric-label">{label}</div>'
        f'<div class="metric-value">{value}</div></div>',
        unsafe_allow_html=True,
    )

left, right = st.columns([1.55, 1])
with left:
    st.markdown("### Difference-colored 3D overlay")
    st.caption(
        "Blue → red encodes 0 → ≥5 Å residual after the one global fit; yellow spokes mark "
        "residues above the adaptive hotspot threshold."
    )
    st.plotly_chart(analysis["overlay"], width="stretch", config={"displaylogo": False})
with right:
    st.markdown("### Chain mapping")
    chain_rows = []
    for item in result.chain_alignments:
        chain_rows.append(
            {
                "Mobile": item.mobile_chain,
                "Reference": item.reference_chain,
                "TM-ref": round(item.tm_reference, 3),
                "RMSD Å": round(item.rmsd, 2),
                "Aligned": item.aligned_length,
                "Coverage": f"{100 * item.coverage_reference:.1f}%",
            }
        )
    st.dataframe(pd.DataFrame(chain_rows), hide_index=True, width="stretch")
    if result.unmatched_mobile_chains or result.unmatched_reference_chains:
        st.warning(
            "Unmatched chains — mobile: "
            f"{', '.join(result.unmatched_mobile_chains) or 'none'}; reference: "
            f"{', '.join(result.unmatched_reference_chains) or 'none'}."
        )
    for warning in result.warnings:
        st.caption(f"⚠ {warning}")

st.markdown("### Difference localization")
difference = result.difference_summary
if difference:
    difference_columns = st.columns(4)
    difference_metrics = [
        ("Median residual", f"{difference.median_displacement:.2f} Å"),
        ("90th percentile", f"{difference.p90_displacement:.2f} Å"),
        ("Largest residual", f"{difference.max_displacement:.2f} Å"),
        (
            f"Residues ≥ {difference.hotspot_threshold:.2f} Å",
            f"{difference.residues_above_threshold} ({100 * difference.fraction_above_threshold:.1f}%)",
        ),
    ]
    for column, (label, value) in zip(difference_columns, difference_metrics, strict=True):
        column.metric(label, value)
    st.markdown(
        '<div class="difference-scale"></div><div class="difference-scale-labels">'
        '<span>near · 0 Å</span><span>moderate · 2 Å</span><span>far · ≥5 Å</span></div>',
        unsafe_allow_html=True,
    )

if result.difference_regions:
    hotspot_rows = [
        {
            "Rank": region.rank,
            "Chain map": f"{region.mobile_chain} → {region.reference_chain}",
            "Mobile region": f"{region.mobile_start}–{region.mobile_end}",
            "Reference region": f"{region.reference_start}–{region.reference_end}",
            "Mean Å": round(region.mean_displacement, 2),
            "Peak Å": round(region.peak_displacement, 2),
            "High residues": f"{region.residues_above_threshold}/{region.residue_count}",
        }
        for region in result.difference_regions
    ]
    st.dataframe(pd.DataFrame(hotspot_rows), hide_index=True, width="stretch")
else:
    st.success("No matched residues exceeded the adaptive difference threshold.")
st.caption(
    "Regions group nearby high-residual residues and add two aligned residues of context. "
    "They localize global-fit disagreement; they do not by themselves establish flexibility."
)

st.markdown("### Per-residue displacement")
st.caption("Every value uses the same global assembly transform; chains are not independently refit.")
st.plotly_chart(analysis["displacement"], width="stretch", config={"displaylogo": False})

st.markdown("### Gained and lost contacts")
contact = result.contact_summary
if contact:
    contact_columns = st.columns(4)
    contact_metrics = [
        ("Gained in mobile", str(contact.gained_contacts)),
        ("Lost from mobile", str(contact.lost_contacts)),
        ("Inter-chain gained", str(contact.inter_chain_gained_contacts)),
        ("Inter-chain lost", str(contact.inter_chain_lost_contacts)),
    ]
    for column, (label, value) in zip(contact_columns, contact_metrics, strict=True):
        column.metric(label, value)
st.caption(
    "Representative-atom contact ≤8 Å; a 1 Å transition margin suppresses threshold jitter. "
    "Same-chain neighbors within two sequence positions are excluded. Larger symbols are inter-chain."
)
contact_chart, contact_table = st.columns([1.15, 1])
with contact_chart:
    st.plotly_chart(analysis["contacts"], width="stretch", config={"displaylogo": False})
with contact_table:
    contact_rows = [
        {
            "Change": row.change,
            "Mobile pair": (
                f"{row.mobile_chain_1}:{row.mobile_residue_1} — "
                f"{row.mobile_chain_2}:{row.mobile_residue_2}"
            ),
            "Reference pair": (
                f"{row.reference_chain_1}:{row.reference_residue_1} — "
                f"{row.reference_chain_2}:{row.reference_residue_2}"
            ),
            "Mobile Å": round(row.mobile_distance, 2),
            "Reference Å": round(row.reference_distance, 2),
            "Inter-chain": row.inter_chain,
        }
        for row in result.contact_changes[:30]
    ]
    if contact_rows:
        st.dataframe(pd.DataFrame(contact_rows), hide_index=True, width="stretch", height=430)
    else:
        st.success("No stable gained or lost contacts were detected.")

st.markdown("### Rigid-body motion vs internal deformation")
st.caption(
    "Each mapped chain is diagnostically refitted after the complex-level fit. The improvement is "
    "the chain rigid-body component; residual after that local fit is internal deformation."
)
motion_chart, motion_table = st.columns([1.05, 1])
with motion_chart:
    st.plotly_chart(analysis["motion"], width="stretch", config={"displaylogo": False})
with motion_table:
    motion_rows = [
        {
            "Chain map": f"{row.mobile_chain} → {row.reference_chain}",
            "Global Å": round(row.global_rmsd, 2),
            "Rigid-body Å": round(row.rigid_body_rmsd, 2),
            "Internal Å": round(row.internal_rmsd, 2),
            "Centroid shift Å": round(row.centroid_shift, 2),
            "Correction °": round(row.correction_rotation_degrees, 2),
            "Pattern": row.classification,
        }
        for row in result.motion_decomposition
    ]
    st.dataframe(pd.DataFrame(motion_rows), hide_index=True, width="stretch")
st.caption(
    "The two RMSD components combine by root-sum-square, not ordinary addition. "
    "The local fit is diagnostic only and never replaces the global residual map."
)

st.markdown("### 3D spatial difference patches")
st.caption(
    "High-residual residues within 10 Å in the reference frame are connected into patches, "
    "even when they are distant in sequence or lie on different chains. P labels appear in the 3D view."
)
patch_rows = [
    {
        "Patch": f"P{patch.rank}",
        "Residues": patch.member_count,
        "Mobile chains": patch.mobile_chains,
        "Reference chains": patch.reference_chains,
        "Mean Å": round(patch.mean_displacement, 2),
        "Peak Å": round(patch.peak_displacement, 2),
        "Radius Å": round(patch.spatial_radius, 2),
        "Mobile members": patch.mobile_residues,
    }
    for patch in result.spatial_patches
]
if patch_rows:
    st.dataframe(pd.DataFrame(patch_rows), hide_index=True, width="stretch")
else:
    st.success("No spatial patches exceeded the adaptive difference threshold.")

st.markdown("### Alignment and unmatched-chain track")
st.markdown(gap_track_html(result), unsafe_allow_html=True)

card_column, download_column = st.columns([1.35, 1])
with card_column:
    st.markdown("### Publication card")
    card_svg = analysis["files"]["card_svg"].decode("utf-8")
    st.markdown(f'<div class="publication-card">{card_svg}</div>', unsafe_allow_html=True)
with download_column:
    st.markdown("### Download result bundle")
    download_specs = [
        ("Residue mapping CSV", "residue_csv", "residue_mapping.csv", "text/csv"),
        ("Chain mapping CSV", "chain_csv", "chain_mapping.csv", "text/csv"),
        ("Difference hotspots CSV", "hotspots_csv", "difference_hotspots.csv", "text/csv"),
        ("Contact changes CSV", "contact_changes_csv", "contact_changes.csv", "text/csv"),
        ("Contact summary JSON", "contact_summary_json", "contact_summary.json", "application/json"),
        ("Motion decomposition CSV", "motion_csv", "motion_decomposition.csv", "text/csv"),
        ("Spatial patches CSV", "spatial_patches_csv", "spatial_patches.csv", "text/csv"),
        (
            "Difference summary JSON",
            "difference_summary_json",
            "difference_summary.json",
            "application/json",
        ),
        ("Publication card SVG", "card_svg", "publication_card.svg", "image/svg+xml"),
        ("Displacement plot SVG", "displacement_svg", "per_residue_displacement.svg", "image/svg+xml"),
        ("Reference PDB", "reference_pdb", "reference.pdb", "chemical/x-pdb"),
        ("Aligned mobile PDB", "aligned_pdb", "mobile_aligned_to_reference.pdb", "chemical/x-pdb"),
        (
            "Difference-colored PDB",
            "difference_pdb",
            "mobile_aligned_difference.pdb",
            "chemical/x-pdb",
        ),
        ("PyMOL difference view", "pymol_script", "view_difference.pml", "text/plain"),
        ("Complete result JSON", "result_json", "result.json", "application/json"),
    ]
    for label, key, filename, mime in download_specs:
        st.download_button(
            label,
            analysis["files"][key],
            file_name=filename,
            mime=mime,
            width="stretch",
        )

with st.expander("Reproducibility and raw US-align output"):
    st.markdown(
        f'<div class="provenance">Engine: US-align {result.engine_version}<br>'
        f'Command: {" ".join(result.command)}<br>'
        f'Mobile SHA-256: {result.mobile_sha256}<br>'
        f'Reference SHA-256: {result.reference_sha256}</div>',
        unsafe_allow_html=True,
    )
    st.code(result.raw_stdout, language="text")
