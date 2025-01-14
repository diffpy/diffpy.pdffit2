.. _examples:

Examples
########

Welcome! This guide offers several examples to help you effectively utilize this package.

Files needed:

    1. :download:`Ni-xray.gr <examples/Ni-xray.gr>` - experimental X-ray PDF data
    2. :download:`Ni.stru <examples/Ni.stru>` - Ni f.c.c. structure in PDFfit format

======================================
Example 1: Calculate PDF of FCC nickel
======================================

The first example shows how to calculate the PDF for FCC nickel and save the resulting data
to a file and plot it using matplotlib.

1. Import the ``PdfFit`` class from the ``diffpy.pdffit2`` module

.. code-block:: python

    from diffpy.pdffit2 import PdfFit

2. Create a PDF calculator object and assign it to the variable ``P``.
Make sure the ``Ni.stru`` file is in the same directory as the script (or specify
the relative path to where it resides)
and you are currently in this directory. Then use ``read_struct`` to read the structure
file, ``alloc`` to configure the PDF
calculation, and then use ``calc`` to carry out the the calculation.

.. code-block:: python

    # create new PDF calculator object
    P = PdfFit()

    # load structure file in PDFFIT or DISCUS format
    P.read_struct("Ni.stru")

    radiation_type = "X"  # x-rays
    qmax = 30.0  # Q-cutoff used in PDF calculation in 1/A
    qdamp = 0.01  # instrument Q-resolution factor, responsible for PDF decay
    rmin = 0.01  # minimum r-value
    rmax = 30.0  # maximum r-value
    npts = 3000  # number of points in the r-grid

    # allocate and configure PDF calculation
    P.alloc(radiation_type, qmax, qdamp, rmin, rmax, npts)
    P.calc()

3. Save the calculated PDF.

.. code-block:: python

    P.save_pdf(1, "Ni_calculation.cgr")

4. We can also plot it using matplotlib

.. code-block:: python

    import matplotlib.pyplot as plt

    # obtain list of r-points and corresponding G values
    r = P.getR()
    G = P.getpdf_fit()

    plt.plot(r, G)
    plt.xlabel("r (Å)")
    plt.ylabel("G (Å$^{-2}$)")
    plt.title("x-ray PDF of nickel simulated at Qmax = %g" % qmax)

    # display plot window, this must be the last command in the script
    plt.show()

The scripts used in this example can be
downloaded :download:`here <examples/Ni_calculation.py>`.

=======================================
Example 2: Performing simple refinement
=======================================

The second example shows how to perform a simple refinement of the Ni structure to
the experimental x-ray PDF. The example uses the same data files as the first example.

1. Import the PdfFit class from the diffpy.pdffit2 module

.. code-block:: python

    from diffpy.pdffit2 import PdfFit

2. Create a PDF calculator object and assign it to the variable ``pf``.
Load the experimental x-ray PDF data using ``read_data`` and also load
the nickel structure file.

.. code-block:: python

    # Create new PDF calculator object.
    pf = PdfFit()

    # Load experimental x-ray PDF data
    qmax = 30.0  # Q-cutoff used in PDF calculation in 1/A
    qdamp = 0.01  # instrument Q-resolution factor, responsible for PDF decay
    pf.read_data("Ni-xray.gr", "X", qmax, qdamp)

    # Load nickel structure, must be in PDFFIT or DISCUS format
    pf.read_struct("Ni.stru")

3. Configure the refinement, assigning structural parameters to variables.  For more
information on how to do this correctly, please read the PDFgui documentation.
Set initial values for the variables using ``setpar``.
Finally, you can configure the range over which to refine (``pdfrange``) and
run the refinement (``refine``).

.. code-block:: python

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

    # Refine all parameters
    pf.pdfrange(1, 1.5, 19.99)
    pf.refine()

4. Save the refined result. ``save_struct`` saves the new, refined, structure to a
``.stru`` format file.  ``save_res`` saves the outputs in a structured text file.

.. code-block:: python

    pf.save_pdf(1, "Ni_refinement.fgr")
    pf.save_struct(1, "Ni_refinement.rstr")
    pf.save_res("Ni_refinement.res")

5. We can also plot it using matplotlib

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy

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

The scripts can be downloaded from :download:`here <examples/Ni_refinement.py>`.
