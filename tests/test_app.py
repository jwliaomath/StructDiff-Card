from pathlib import Path

from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]


def test_bundled_example_renders_without_exception() -> None:
    app = AppTest.from_file(str(ROOT / "app.py"), default_timeout=30).run()
    assert not app.exception

    example_button = next(
        button for button in app.button if button.label.startswith("Load 2HHB")
    )
    app = example_button.click().run(timeout=30)

    assert not app.exception
    assert any(markdown.value == "### Publication card" for markdown in app.markdown)
    assert any(markdown.value == "### Gained and lost contacts" for markdown in app.markdown)
    assert any(
        markdown.value == "### Rigid-body motion vs internal deformation"
        for markdown in app.markdown
    )
    assert any(
        markdown.value == "### 3D spatial difference patches" for markdown in app.markdown
    )
