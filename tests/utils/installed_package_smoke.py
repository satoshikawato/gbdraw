"""Run with a fresh venv's Python from a directory containing copied inputs."""

from pathlib import Path
import importlib.metadata
import json
import subprocess
import sys
import sysconfig
import xml.etree.ElementTree as ET

import gbdraw


def main() -> None:
    package_path = Path(gbdraw.__file__).resolve()
    assert package_path.is_relative_to(Path(sysconfig.get_path("purelib")).resolve())
    assert sys.prefix != sys.base_prefix
    assert gbdraw.__version__ == importlib.metadata.version("gbdraw")
    # App-shell startup fetches this JSON before the first Generate.
    palettes = json.loads((package_path.parent / "web/gallery/palettes/palettes.json").read_text())
    assert palettes["palettes"]["default"]["CDS"]
    cli = Path(sysconfig.get_path("scripts")) / (
        "gbdraw.exe" if sys.platform == "win32" else "gbdraw"
    )

    def run(*args: str) -> str:
        result = subprocess.run([str(cli), *args], check=True, capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)
        return result.stdout

    assert gbdraw.__version__ in run("--version")
    assert "circular" in run("--help")
    for mode in ("circular", "linear"):
        assert "--gbk" in run(mode, "--help")
        run(mode, "--gbk", "HmmtDNA.gbk", "-o", mode,
            "--session_output", f"{mode}.session.json")
        run(mode, "--session", f"{mode}.session.json", "-o", f"{mode}-replay")
        original = ET.parse(f"{mode}.svg").getroot()
        replay = ET.parse(f"{mode}-replay.svg").getroot()
        assert original.tag == "{http://www.w3.org/2000/svg}svg"
        assert original.findall(".//{http://www.w3.org/2000/svg}path")
        assert original.findall(".//{http://www.w3.org/2000/svg}text")
        assert ET.tostring(original) == ET.tostring(replay)

    records = gbdraw.read_genbank("HmmtDNA.gbk")
    for mode, diagram in (
        ("circular", gbdraw.draw_circular(records[0])),
        ("linear", gbdraw.draw_linear(records)),
    ):
        diagram.save(f"api-{mode}.svg")
        assert ET.parse(f"api-{mode}.svg").findall(".//{http://www.w3.org/2000/svg}path")

    print(json.dumps({
        "python": sys.version, "prefix": sys.prefix,
        "package_path": str(package_path), "version": gbdraw.__version__,
        "sys_path": sys.path, "result": "PASS",
    }, indent=2))


if __name__ == "__main__":
    main()
