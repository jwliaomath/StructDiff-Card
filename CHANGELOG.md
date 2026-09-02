# Changelog

## 0.3.0 — 2026-09-02

- Add stable gained/lost representative-atom contact analysis and contact map.
- Add per-chain rigid-body versus internal-deformation diagnostics.
- Add 3D spatial clustering of high-residual residues across sequence and chain boundaries.
- Reduce local 3D residue markers and repair sidebar contrast for uploads and buttons.
- Export contact, motion-decomposition, and spatial-patch result tables.

## 0.2.0 — 2026-09-02

- Difference-colored 3D overlay with a fixed 0–5 Å scale and high-residual vectors.
- Adaptive coherent-region detection using the larger of 2 Å and the 90th percentile.
- Seven-residue smoothed displacement trace plus numbered region peaks.
- Difference summary JSON and coherent-region CSV exports.
- Aligned PDB with residual displacement stored in B-factor and a ready-to-run PyMOL view.

## 0.1.0 — 2026-09-02

- Multi-chain US-align dispatcher with automatic or manual chain mapping.
- One global assembly transform for per-residue displacement.
- Explicit unmatched-chain and sequence-gap tracks.
- Streamlit interface, CLI, Python API, Docker image, and reproducible exports.
- Static GitHub Pages gallery with a precomputed 2HHB/1HHO result and 15-second GIF.