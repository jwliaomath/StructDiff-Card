"""Exercise upload staging without importing Streamlit's top-level UI."""
import ast
import tempfile
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "mobile_name, reference_name",
    [
        ("10.pdb", "10.pdb"),
        ("model.pdb", "MODEL.PDB"),
        ("mobile.pdb", "reference.cif"),
    ],
)
def test_run_analysis_keeps_both_input_payloads(mobile_name, reference_name):
    source = Path(__file__).resolve().parents[1] / "app.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "run_analysis")
    module = ast.Module(body=[function], type_ignores=[])

    class StopAfterInputCheck(Exception):
        pass

    def inspect_inputs(mobile, reference, **kwargs):
        assert mobile != reference
        assert mobile.parent != reference.parent
        assert mobile.name == mobile_name
        assert reference.name == reference_name
        assert mobile.read_bytes() == b"mobile coordinates"
        assert reference.read_bytes() == b"reference coordinates"
        raise StopAfterInputCheck

    namespace = {
        "Path": Path,
        "tempfile": tempfile,
        "compare_structures": inspect_inputs,
        "parse_manual_mapping": lambda value: None,
    }
    exec(compile(module, str(source), "exec"), namespace)  # noqa: S102
    with pytest.raises(StopAfterInputCheck):
        namespace["run_analysis"](
            b"mobile coordinates",
            mobile_name,
            b"reference coordinates",
            reference_name,
            "asymmetric",
            "",
        )
