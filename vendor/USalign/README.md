# Vendored US-align source

- Upstream: https://zhanggroup.org/US-align/bin/module/USalign.cpp
- Official repository: https://github.com/pylelab/USalign
- Version printed by this source: `20260527`
- Upstream changelog embedded in the snapshot: through `2026-08-26`
- Downloaded: 2026-09-02
- SHA-256: `5e05ddcd66d7ca7553954b7df94e938fcfbf1a5d94a81edf5b04e80f8d8913f5`

The complete upstream notice, references, and warranty disclaimer are retained
at the head of `USalign.cpp`. StructDiff Card does not modify this source. The
Dockerfile verifies the checksum before compiling it.

To update US-align, replace the source, update the checksum in this file,
`Dockerfile`, and `.github/workflows/tests.yml`, then regenerate the example and
review parser fixtures against the new output.
