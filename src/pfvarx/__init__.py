"""PfVarX package.

``__version__`` below is the single authoritative pipeline release version. The package
metadata reads it through ``[tool.hatch.version]`` in pyproject.toml, and every report,
VCF header, snapshot and source archive takes it from here. Change it only here; the test
``tests/unit/test_version.py`` fails if any current-release reference disagrees.
"""

__version__ = "0.0.0"
