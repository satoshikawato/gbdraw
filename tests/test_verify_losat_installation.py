import ctypes
import errno
from types import SimpleNamespace

import pytest

from tools import verify_losat_installation as acceptance


@pytest.mark.parametrize("value,error,expected_mode,exception", [
    (1, 0, "rosetta", None),
    (0, 0, "native", None),
    (0, errno.ENOENT, "native", None),
    (0, 0, "rosetta", ValueError),
    (1, 0, "native", ValueError),
    (0, errno.EPERM, "native", OSError),
    (2, 0, "native", ValueError),
])
def test_rosetta_requires_current_process_evidence(monkeypatch, value, error, expected_mode, exception):
    monkeypatch.setattr(acceptance.platform, "system", lambda: "Darwin")
    # x86_64 alone cannot distinguish native Intel from translated Apple Silicon.
    monkeypatch.setattr(acceptance.platform, "machine", lambda: "x86_64")

    def sysctl(name, output, size, new_value, new_size):
        assert name == b"sysctl.proc_translated"
        ctypes.cast(output, ctypes.POINTER(ctypes.c_int))[0] = value
        ctypes.set_errno(error)
        return -1 if error else 0

    monkeypatch.setattr(acceptance.ctypes, "CDLL", lambda *a, **kw: SimpleNamespace(sysctlbyname=sysctl))
    if exception:
        with pytest.raises(exception):
            acceptance.runtime_identity(expected_mode)
    else:
        assert acceptance.runtime_identity(expected_mode)["execution_mode"] == expected_mode


def test_non_macos_cannot_claim_rosetta(monkeypatch):
    monkeypatch.setattr(acceptance.platform, "system", lambda: "Linux")
    with pytest.raises(ValueError, match="Expected rosetta"):
        acceptance.runtime_identity("rosetta")


def test_macos_metadata_uses_read_only_native_xattr(monkeypatch, tmp_path):
    monkeypatch.setattr(acceptance.platform, "system", lambda: "Darwin")
    binary = tmp_path / "LOSAT"
    calls = []

    def output(argv, **kwargs):
        calls.append(argv)
        if argv == ["/usr/bin/xattr", str(binary)]:
            return "com.apple.quarantine\n"
        assert argv == ["/usr/bin/xattr", "-p", "-x", "com.apple.quarantine", str(binary)]
        return "30 30 38 31\n"

    monkeypatch.setattr(acceptance.subprocess, "check_output", output)
    monkeypatch.setattr(acceptance.subprocess, "run", lambda *a, **kw:
                        SimpleNamespace(returncode=1, stdout="", stderr="not signed"))
    metadata = acceptance.macos_binary_metadata(binary)
    assert metadata["quarantine_present"]
    assert metadata["extended_attributes_hex"] == {"com.apple.quarantine": "30303831"}
    assert len(calls) == 2
