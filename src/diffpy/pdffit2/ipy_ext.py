#!/usr/bin/env python

"""This module defines functions within IPython session to simulate the
old pdffit2 interactive session.

Usage: %load_ext diffpy.pdffit2.ipy_ext
"""

from diffpy.utils._deprecator import build_deprecation_message, deprecated

base = "diffpy.pdffit2.ipy_ext.EasyPDFPlotting"
removal_version = "2.0.0"
showRw_deprecation_msg = build_deprecation_message(
    base,
    "showRw",
    "show_rw",
    removal_version,
)

# extension begin


def load_ipython_extension(ipython):
    from diffpy.pdffit2 import PdfFit

    pf = PdfFit()
    pdf = EasyPDFPlotting(pf)
    print("        Type  help(pdffit)  or  help(topic)  for information.\n")
    ns = dict(pdffit=PdfFit, pdf=pdf)
    pf._export_all(ns)
    ipython.user_ns.update(ns)
    return


class EasyPDFPlotting(object):
    """Convenience functions for accessing and plotting PDFfit2 data."""

    def __init__(self, pdffit_instance):
        self._pdffit = pdffit_instance
        return

    @property
    def r(self):
        "R-grid for PDF simulation."
        return self._asarray(self._pdffit.get_r(), dtype=float)

    @property
    def Gobs(self):
        "Observed PDF data."
        return self._asarray(self._pdffit.getpdf_obs(), dtype=float)

    @property
    def Gcalc(self):
        "Calculated PDF data."
        return self._asarray(self._pdffit.getpdf_fit())

    @property
    def Gdiff(self):
        "Difference between the observed and simulated PDF."
        return self.Gobs - self.Gcalc

    def showfit(self, offset=None):
        """Plot observed and simulated PDFs and the difference curve.

        offset   -- offset for the difference curve.

        No return value.
        """
        from math import floor

        from matplotlib.pyplot import gca

        cr = self.r
        cGobs = self.Gobs
        cGcalc = self.Gcalc
        cGdiff = self.Gdiff
        if offset is None:
            offset = floor(min([min(cGobs), min(cGcalc)]) - max(cGdiff))
        ax = gca()
        ax.plot(cr, cGobs, "r.", cr, cGcalc, "b-", cr, cGdiff + offset, "g-")
        xl = ax.xaxis.get_label().get_text()
        yl = ax.yaxis.get_label().get_text()
        if xl == "":
            ax.set_xlabel("r (A)")
        if yl == "":
            ax.set_ylabel("G (A**-2)")
        return

    def show_rw(self):
        """Plot cumulative Rw."""
        from matplotlib.pyplot import gca

        cRw = self._asarray(self._pdffit.getcrw())
        ax = gca()
        ax.plot(self.r, cRw)
        ax.set_title("Cumulative Rw = %.4f" % cRw[-1])
        ax.set_xlabel("r")
        ax.set_ylabel("Rw")
        return

    @deprecated(showRw_deprecation_msg)
    def showRw(self):
        """This function has been deprecated and will be removed in
        version 2.0.0.

        Please use diffpy.pdffit2.ipy_ext.EasyPDFPlotting.show_rw
        instead.
        """
        return self.show_rw()

    @staticmethod
    def _asarray(x, dtype=None):
        import numpy

        return numpy.asarray(x, dtype=dtype)


# End of class EasyPDFPlotting
