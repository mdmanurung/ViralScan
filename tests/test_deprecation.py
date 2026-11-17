"""DEPRECATE-01: a deprecated interface names its replacement and the first release that may remove it."""

import pytest

from viralscan.cli.parsers import parse_args
from viralscan.sensitivity import fragment_capture, fragment_capture_exact
from viralscan.utils import DEPRECATION_REMOVAL, deprecation_message, warn_deprecated


def test_message_names_replacement_and_removal_release() -> None:
    text = deprecation_message("--old", "--new")
    assert "--old is deprecated" in text and "--new" in text and DEPRECATION_REMOVAL in text


def test_warn_deprecated_is_a_future_warning() -> None:
    with pytest.warns(FutureWarning, match="use b"):
        warn_deprecated("a", "b")


def test_fragment_capture_warns_and_still_delegates() -> None:
    with pytest.warns(FutureWarning, match="fragment_capture_exact"):
        old = fragment_capture(0.1, read_length=90)
    assert old == fragment_capture_exact(0.1, read_length=90)


def test_bare_form_warns_with_the_replacement_and_removal(capsys) -> None:
    import contextlib

    with contextlib.suppress(SystemExit):
        parse_args(["-o", "out", "--r1", "a_R1.fq", "--r2", "a_R2.fq"])
    err = capsys.readouterr().err
    assert "without a subcommand is deprecated" in err
    assert "viralscan run" in err and DEPRECATION_REMOVAL in err
