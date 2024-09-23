#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Calculate PDF of FCC nickel. Save data to Ni_calculation.cgr and plot it using matplotlib.
"""

import matplotlib.pyplot as plt

from diffpy.pdffit2 import PdfFit

# create new PDF calculator object
P = PdfFit()

# Load data ------------------------------------------------------------------

# load structure file in PDFFIT or DISCUS format
P.read_struct("Ni.stru")

# Configure calculation ------------------------------------------------------

radiation_type = "X"  # x-rays
qmax = 30.0  # Q-cutoff used in PDF calculation in 1/A
qdamp = 0.01  # instrument Q-resolution factor, responsible for PDF decay
rmin = 0.01  # minimum r-value
rmax = 30.0  # maximum r-value
npts = 3000  # number of points in the r-grid

# allocate and configure PDF calculation
P.alloc(radiation_type, qmax, qdamp, rmin, rmax, npts)

# Calculate -------------------------------------------------------------------

P.calc()

# Save results ---------------------------------------------------------------

P.save_pdf(1, "Ni_calculation.cgr")

# Plot results ---------------------------------------------------------------

# obtain list of r-points and corresponding G values
r = P.getR()
G = P.getpdf_fit()

plt.plot(r, G)
plt.xlabel("r (Å)")
plt.ylabel("G (Å$^{-2}$)")
plt.title("x-ray PDF of nickel simulated at Qmax = %g" % qmax)

# display plot window, this must be the last command in the script
plt.show()
