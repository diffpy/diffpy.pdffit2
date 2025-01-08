"""Unit tests for __version__.py."""

import diffpy.pdffit2


def test_package_version():
    """Ensure the package version is defined and not set to the initial
    placeholder."""
    assert hasattr(diffpy.pdffit2, "__version__")
    assert diffpy.pdffit2.__version__ != "0.0.0"
