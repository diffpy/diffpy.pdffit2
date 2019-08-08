# Release notes

Notable differences from version 1.1.

## Unreleased – Version FIXME

### Added

- Support for Python 3.7, 3.6, 3.5 in addition to 2.7.
- Makefile targets `module` and `test`.
- Support installation from a tagged git-archive tarball.

### Changed

- Build Anaconda package with Anaconda C++ compiler.
- Allow language standard c++11.
- Update to diffpy.structure 3.0.

### Deprecated

- Variable `__gitsha__` in the `version` module which was renamed
  to `__git_commit__`.

### Removed

- Support for Python 2.6.

### Fixed

- Windows build with MSVC 9 and MSVC 14.
- Invalid escape sequences in string values.
- Open files within the `with` context so they get closed when done.
