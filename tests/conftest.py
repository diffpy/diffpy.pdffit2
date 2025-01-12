import io
import json
from pathlib import Path

import pytest

import diffpy.pdffit2
import diffpy.pdffit2.output  # assuming this is the correct import path


@pytest.fixture
def user_filesystem(tmp_path):
    base_dir = Path(tmp_path)
    home_dir = base_dir / "home_dir"
    home_dir.mkdir(parents=True, exist_ok=True)
    cwd_dir = base_dir / "cwd_dir"
    cwd_dir.mkdir(parents=True, exist_ok=True)

    home_config_data = {"username": "home_username", "email": "home@email.com"}
    with open(home_dir / "diffpyconfig.json", "w") as f:
        json.dump(home_config_data, f)

    yield tmp_path


@pytest.fixture
def datafile():
    """Fixture to dynamically load any test file."""

    def _load(filename):
        return "tests/testdata/" + filename

    return _load


@pytest.fixture
def capture_output():
    """Capture output from pdffit2 engine produced in function call."""

    def _capture(f, *args, **kwargs):
        savestdout = diffpy.pdffit2.output.stdout
        fp = io.StringIO()
        diffpy.pdffit2.redirect_stdout(fp)
        try:
            f(*args, **kwargs)
        finally:
            diffpy.pdffit2.redirect_stdout(savestdout)
            diffpy.pdffit2.output.restore_stdout()
        return fp.getvalue()

    return _capture
