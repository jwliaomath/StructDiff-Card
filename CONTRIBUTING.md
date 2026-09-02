# Contributing

Small, reviewable contributions are welcome. Please open an issue before a
large feature so scientific scope and test cases can be agreed first.

## Development

```bash
uv sync --extra dev
export USALIGN_BIN=/path/to/USalign
uv run pytest
```

The optional integration test uses the bundled public 2HHB/1HHO example and is
enabled only when `USALIGN_BIN` is set. Frontend changes must also pass:

```bash
cd site
pnpm install --frozen-lockfile
pnpm build
```

Please add a regression fixture for changes to US-align output parsing. Do not
commit private or unpublished coordinate files.