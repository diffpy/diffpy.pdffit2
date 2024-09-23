#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Perform simple refinement of Ni structure to the experimental x-ray PDF.
Save fitted curve, refined structure and results summary.
"""

import matplotlib.pyplot as plt
import numpy

from diffpy.pdffit2 import PdfFit

# Create new PDF calculator object.
pf = PdfFit()

# Load data ------------------------------------------------------------------

# Load experimental x-ray PDF data
qmax = 30.0  # Q-cutoff used in PDF calculation in 1/A
qdamp = 0.01  # instrument Q-resolution factor, responsible for PDF decay
pf.read_data("Ni-xray.gr", "X", qmax, qdamp)

# Load nickel structure, must be in PDFFIT or DISCUS format
pf.read_struct("Ni.stru")

# Configure Refinement -------------------------------------------------------

# Refine lattice parameters a, b, c.
# Make them all equal to parameter @1.
pf.constrain(pf.lat(1), "@1")
pf.constrain(pf.lat(2), "@1")
pf.constrain(pf.lat(3), "@1")
# set initial value of parameter @1
pf.setpar(1, pf.lat(1))

# Refine phase scale factor.  Right side can have formulas.
pf.constrain("pscale", "@20 * 2")
pf.setpar(20, pf.getvar(pf.pscale) / 2.0)

# Refine PDF damping due to instrument Q-resolution.
# Left side can be also passed as a reference to PdfFit object
pf.constrain(pf.qdamp, "@21")
pf.setpar(21, 0.03)

# Refine sharpening factor for correlated motion of close atoms.
pf.constrain(pf.delta2, 22)
pf.setpar(22, 0.0003)

# Set all temperature factors isotropic and equal to @4
for idx in range(1, 5):
    pf.constrain(pf.u11(idx), "@4")
    pf.constrain(pf.u22(idx), "@4")
    pf.constrain(pf.u33(idx), "@4")
pf.setpar(4, pf.u11(1))

# Refine ---------------------------------------------------------------------

pf.pdfrange(1, 1.5, 19.99)
pf.refine()

# Save results ---------------------------------------------------------------

pf.save_pdf(1, "Ni_refinement.fgr")
pf.save_struct(1, "Ni_refinement.rstr")
pf.save_res("Ni_refinement.res")

# Plot results ---------------------------------------------------------------

# obtain data from PdfFit calculator object
r = pf.getR()
Gobs = pf.getpdf_obs()
Gfit = pf.getpdf_fit()

# calculate difference curve
Gdiff = numpy.array(Gobs) - numpy.array(Gfit)
Gdiff_baseline = -10

plt.plot(r, Gobs, "ko")
plt.plot(r, Gfit, "b-")
plt.plot(r, Gdiff + Gdiff_baseline, "r-")

plt.xlabel("r (Å)")
plt.ylabel("G (Å$^{-2}$)")
plt.title("Fit of nickel to x-ray experimental PDF")

# display plot window, this must be the last command in the script
plt.show()
