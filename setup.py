from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

from setuptools import find_packages, setup


def _load_build_support_module():
    # Load the helper by path so PEP 517 isolated builds do not need `gbdraw` importable first.
    build_support_path = Path(__file__).resolve().parent / "gbdraw" / "_build_support.py"
    spec = spec_from_file_location("gbdraw_build_support", build_support_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load build support module from {build_support_path}")

    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_build_support = _load_build_support_module()
get_excluded_package_data_patterns = _build_support.get_excluded_package_data_patterns
get_package_data_patterns = _build_support.get_package_data_patterns
is_browser_wheel_build = _build_support.is_browser_wheel_build


INCLUDE_BROWSER_WHEEL = not is_browser_wheel_build()


setup(
    packages=find_packages(include=["gbdraw*"]),
    package_data={"gbdraw": get_package_data_patterns(include_browser_wheel=INCLUDE_BROWSER_WHEEL)},
    exclude_package_data={"gbdraw": get_excluded_package_data_patterns(include_browser_wheel=INCLUDE_BROWSER_WHEEL)},
    include_package_data=False,
)
