########################################################################
#
# pdffit2           by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

"""Definition of __version__ for pdffit2 package.
"""

__id__ = "$Id$"

try:
    from diffpy.version import __svnrevision__, __date__
except ImportError:
    __svnrevision__ = "?"
    __date__ = "?"

__version__ = "2.0." + __svnrevision__

# End of file
