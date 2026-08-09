#!/usr/bin/env python

"""Tests for deprecated public APIs."""

from unittest.mock import Mock

import pytest

from diffpy.pdffit2.ipy_ext import EasyPDFPlotting
from diffpy.pdffit2.pdffit import PdfFit


@pytest.mark.parametrize(
    ("legacy_name", "replacement_name", "args"),
    [
        ("getR", "get_r", ()),
        ("selectAtomType", "select_atom_type", (1, "i", "Ni", True)),
        ("selectAtomIndex", "select_atom_index", (1, "j", 2, False)),
        ("selectAll", "select_all", (1, "i")),
        ("selectNone", "select_none", (1, "j")),
    ],
)
def test_pdffit_camel_case_method_is_deprecated(
    monkeypatch, legacy_name, replacement_name, args
):
    """Camel-case PdfFit methods warn and delegate to snake case."""
    pdf_fit = PdfFit()
    expected = object()
    replacement = Mock(return_value=expected)
    monkeypatch.setattr(pdf_fit, replacement_name, replacement)

    with pytest.warns(DeprecationWarning, match=replacement_name):
        actual = getattr(pdf_fit, legacy_name)(*args)

    assert actual is expected
    replacement.assert_called_once_with(*args)


def test_show_rw_camel_case_method_is_deprecated(monkeypatch):
    """EasyPDFPlotting.showRw warns and delegates to show_rw."""
    plotting = EasyPDFPlotting(Mock())
    expected = object()
    replacement = Mock(return_value=expected)
    monkeypatch.setattr(plotting, "show_rw", replacement)

    with pytest.warns(DeprecationWarning, match="show_rw"):
        actual = plotting.showRw()

    assert actual is expected
    replacement.assert_called_once_with()
