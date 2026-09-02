# Security and privacy

StructDiff Card is designed for local use. The published GitHub Pages site is a
static, precomputed gallery and has no upload endpoint.

- The Python wrapper invokes US-align with an argument list and `shell=False`.
- User-controlled raw US-align flags are not accepted.
- Streamlit analyses use isolated temporary directories; uploaded structures
  are converted to in-memory result objects before the directory is removed.
- The Docker container runs as an unprivileged user and the Compose profile uses
  a read-only root filesystem plus a bounded temporary filesystem.

Please report a vulnerability privately through GitHub Security Advisories if
the repository has them enabled. Do not attach sensitive structure files to a
public issue.