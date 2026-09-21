import ctypes
import errno
from pathlib import Path
import sys
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


def test_extracted_source_identity_requires_declared_commit(tmp_path: Path) -> None:
    source = tmp_path / "LOSAT-0.1.0"
    source.mkdir()
    sha = "6bfb1b09b6cb9451fa771e687c82cbb860e8c779"
    assert acceptance.source_identity(source, sha) == sha
    with pytest.raises(ValueError, match="--losat-source-sha"):
        acceptance.source_identity(source, None)


def test_conda_prepare_uses_running_prefix_without_setup(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    prefix = tmp_path / "conda"
    (prefix / "conda-meta").mkdir(parents=True)
    binary = prefix / "bin" / "losat"
    binary.parent.mkdir()
    binary.write_text("#!/bin/sh\n", encoding="utf-8")
    binary.chmod(0o755)
    original = binary.read_bytes()
    monkeypatch.setattr(acceptance.sys, "prefix", str(prefix))
    monkeypatch.setattr(
        acceptance.setup,
        "setup_losat",
        lambda: (_ for _ in ()).throw(AssertionError("setup must not run")),
    )
    monkeypatch.setattr(
        acceptance.setup,
        "managed_losat",
        lambda: (_ for _ in ()).throw(AssertionError("cache must not be read")),
    )

    runtime = acceptance.prepare_runtime("conda", ["unused-python"])

    assert runtime.source == "conda"
    assert runtime.executable == str(binary.absolute())
    assert binary.read_bytes() == original


def test_conda_prepare_rejects_normal_path_fallback(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    prefix = tmp_path / "conda"
    (prefix / "conda-meta").mkdir(parents=True)
    path_binary = tmp_path / "other-environment" / "losat"
    monkeypatch.setattr(acceptance.sys, "prefix", str(prefix))
    monkeypatch.setattr(acceptance.setup, "managed_losat", lambda: None)
    monkeypatch.setattr(acceptance.protein, "_bundled_losatp_resource", lambda: None)
    monkeypatch.setattr(
        acceptance.protein,
        "_path_executable",
        lambda name: str(path_binary) if name == "losat" else None,
    )

    with pytest.raises(ValueError, match="requires source=conda"):
        acceptance.prepare_runtime("conda", ["unused-python"])


@pytest.mark.parametrize(
    ("system", "machine", "expected"),
    [
        ("Linux", "x86_64", "linux-64"),
        ("Darwin", "x86_64", "osx-64"),
        ("Darwin", "arm64", "osx-arm64"),
    ],
)
def test_conda_subdir_is_limited_to_initial_targets(
    monkeypatch: pytest.MonkeyPatch,
    system: str,
    machine: str,
    expected: str,
) -> None:
    monkeypatch.setattr(acceptance.platform, "system", lambda: system)
    monkeypatch.setattr(acceptance.platform, "machine", lambda: machine)
    assert acceptance.conda_subdir() == expected


def test_python_api_smoke_uses_public_package_root(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    records = [SimpleNamespace(id="one"), SimpleNamespace(id="two")]
    captured = {}

    class FakeDiagram:
        def save(self, path):
            Path(path).write_text("<svg/>", encoding="utf-8")
            return Path(path)

    def fake_draw_linear(actual_records, *, options):
        captured["records"] = actual_records
        captured["options"] = options
        return FakeDiagram()

    monkeypatch.setattr(acceptance.gbdraw, "draw_linear", fake_draw_linear)

    report = acceptance.run_python_api_smoke(records, tmp_path)

    assert captured["records"] is records
    assert captured["options"].comparisons.protein_mode == "pairwise"
    assert captured["options"].comparisons.threads == 1
    assert report["entrypoint"] == "gbdraw.draw_linear"
    assert Path(report["svg"]).read_text(encoding="utf-8") == "<svg/>"


def test_cli_smoke_runs_installed_cli_in_process_for_native_argv_capture(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    records = []
    for index in range(2):
        record = acceptance.SeqRecord(acceptance.Seq("ATG"), id=f"record{index}")
        record.annotations["molecule_type"] = "DNA"
        records.append(record)

    def fake_main() -> None:
        raw_output = Path(sys.argv[sys.argv.index("--protein_blastp_output") + 1])
        output_prefix = Path(sys.argv[sys.argv.index("--output") + 1])
        raw_output.write_text("query\tsubject\n", encoding="utf-8")
        output_prefix.with_suffix(".svg").write_text("<svg/>", encoding="utf-8")

    monkeypatch.setattr(acceptance.gbdraw_cli, "main", fake_main)

    report = acceptance.run_cli_smoke(records, tmp_path)

    assert report["returncode"] == 0
    assert report["comparison_rows"] == 1
    assert report["entrypoint"] == "gbdraw.cli.main"
    assert report["argv"][0] == "gbdraw"
    assert report["replay_command"][:3] == [sys.executable, "-m", "gbdraw.cli"]
