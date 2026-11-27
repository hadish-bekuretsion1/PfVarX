from pathlib import Path

import pytest

from pfvarx.caller_adapter import check_backend
from pfvarx.configuration.models import CallerConfig
from pfvarx.io import sha256


def test_external_backend_requires_matching_executable(tmp_path: Path) -> None:
    executable = tmp_path / "caller adapter"
    executable.write_text("#!/bin/sh\nexit 0\n")
    executable.chmod(0o700)
    config = CallerConfig(backend="external", executable=executable, sha256=sha256(executable))
    check_backend(config)
    executable.write_text("#!/bin/sh\nexit 1\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        check_backend(config)


def test_native_backend_is_default_for_older_config() -> None:
    assert CallerConfig.model_validate({}).backend == "native"
