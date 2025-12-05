=============
Release notes
=============

.. current developments

1.5.2
=====

**Added:**

* Add gsl to conda requirements.

**Changed:**

* Remove mac static GSL linking; always link to shared GSL.
* Change doc to docs for skpkg standard.

**Fixed:**

* Fix misspelled words in source code comments.
* Migrate documentation to `scikit-package 0.1.0` standards, including a mock import for API rendering.
* Add ``libsblas.dll`` to build pypi wheel in windows.


1.5.1
=====

**Fixed:**

* Fixed `SystemError` and `MemoryError` for `redirect_stdout` on Windows with Python 3.13.

**Removed:**

* Removed `restore_stdout` function and wrapper.


1.5.0
=====

**Added:**

* Python 3.11, 3.12 support
* Option to skip printing of introductory information when initializing the `PdfFit` class.
* Added additional runtime linker flags in `CustomBuildExt.run` to embed the `RPATH` flags for the built extensions.
* Support for retrieving GSL configuration from `CONDA_PREFIX`/ `GSL_PATH` on all platforms.
* Separate installation instruction for macOS (Arm64) in READEM
* Added `restore_stdout` function and wrapper.
* Added Python 3.13 support.

**Changed:**

* Changed setup.py to lazy evaluate gsl installation.
* Documentation brought up to date
* Merged the GSL configuration logic in `setup.py`.
* Changed `pytest` `capture_output` fixture. Now automatically restores `sys.stdout`.

**Fixed:**

* remove older conda-recipe files
* moved the tests directory from src to the root using conftest.py.
* fixed a circular import bug during " pip install ." in GitHub CI.
* renamed .py files under tests to snake_case.
* add PyPI packages under pip.txt
* re-cookiecutter to group's package standard
* Fix missing `__date__`, use PyPI release date.
* Fixed `SystemError` when running `pytest` on Windows with Python 3.13.

**Removed:**

* Python <= 3.10 support
* Six dependency and py2 support


1.4.2
=====

**Added:**

* Support for Python 3.11, 3.12

**Changed:**

* No notable functional changes from 1.4.1

1.4.4rc0
========

**Fixed:**

* Code linted to group flake8 standards
* Package structure moved to diffpy standard structure
