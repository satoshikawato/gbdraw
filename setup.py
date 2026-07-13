from pathlib import Path
from runpy import run_path
from types import SimpleNamespace

from setuptools import find_packages, setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.sdist import sdist as _sdist


def _load_build_support_module():
    # Load the helper by path so PEP 517 isolated builds do not need `gbdraw` importable first.
    path = Path(__file__).resolve().parent / "gbdraw" / "_build_support.py"
    return SimpleNamespace(**run_path(str(path)))


_build_support = _load_build_support_module()
get_excluded_package_data_patterns = _build_support.get_excluded_package_data_patterns
get_package_data_patterns = _build_support.get_package_data_patterns
is_browser_wheel_build = _build_support.is_browser_wheel_build


INCLUDE_BROWSER_WHEEL = not is_browser_wheel_build()


class build_py(_build_py):
    def run(self):
        _build_support.refresh_open_source_notices()
        super().run()


class sdist(_sdist):
    def run(self):
        _build_support.refresh_open_source_notices()
        super().run()


setup(
    packages=find_packages(include=["gbdraw*"]),
    package_data={"gbdraw": get_package_data_patterns(include_browser_wheel=INCLUDE_BROWSER_WHEEL)},
    exclude_package_data={"gbdraw": get_excluded_package_data_patterns(include_browser_wheel=INCLUDE_BROWSER_WHEEL)},
    include_package_data=False,
    cmdclass={"build_py": build_py, "sdist": sdist},
)
