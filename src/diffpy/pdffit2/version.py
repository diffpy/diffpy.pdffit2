#!/usr/bin/env python
##############################################################################
#
# (c) 2024 The Trustees of Columbia University in the City of New York.
# All rights reserved.
#
# File coded by: Billinge Group members and community contributors.
#
# See GitHub contributions for a more detailed list of contributors.
# https://github.com/diffpy/diffpy.pdffit2/graphs/contributors
#
# See LICENSE.rst for license information.
#
##############################################################################
"""Definition of __version__."""

import datetime
import json
import urllib.request
from importlib.metadata import version
from pathlib import Path


def get_pypi_release_date(package_name, timeout=5):
    package_file = Path(__file__).resolve()

    try:
        with open(package_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
        for line in reversed(lines):
            if line.startswith("# Release date:"):
                return line.split(":", 1)[1].strip()

        url = f"https://pypi.org/pypi/{package_name}/json"
        with urllib.request.urlopen(url, timeout=timeout) as response:
            data = json.loads(response.read().decode("utf-8"))

        installed_version = version(package_name)
        release_data = data["releases"].get(installed_version, [])
        if not release_data:
            raise ValueError(f"No release data found for version {installed_version}")

        release_date_str = release_data[-1]["upload_time"]
        release_date = datetime.datetime.fromisoformat(release_date_str).date()

        with open(package_file, "a", encoding="utf-8") as f:
            f.write(f"\n# Release date: {release_date}")

    except (ValueError, OSError) as e:
        print(f"Warning: Could not fetch release date from PyPI: {e}")
        release_date = datetime.datetime.fromtimestamp(package_file.stat().st_ctime).isoformat()

    return str(release_date)


__version__ = version("diffpy.pdffit2")
__date__ = get_pypi_release_date("diffpy.pdffit2")

# End of file
