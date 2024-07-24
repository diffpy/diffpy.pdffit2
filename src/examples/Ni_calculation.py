#!/usr/bin/env python

'''Calculate PDF of FCC nickel.  Save data to Ni_calculation.cgr.
'''

from diffpy.pdffit2 import PdfFit

# create new PDF calculator object
P = PdfFit()

# load structure file in PDFFIT or DISCUS format
P.read_struct("Ni.stru")

radiation_type = 'X'  # x-rays
qmax = 30.0  # Q-cutoff used in PDF calculation in 1/A
qdamp = 0.01 # instrument Q-resolution factor, responsible for PDF decay
rmin = 0.01  # minimum r-value
rmax = 30.0  # maximum r-value
npts = 3000  # number of points in the r-grid

# allocate and configure PDF calculation
P.alloc(radiation_type, qmax, qdamp, rmin, rmax, npts)
P.calc()

P.save_pdf(1, "Ni_calculation.cgr")
