from setuptools import find_packages, setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.sdist import sdist as _sdist

from gbdraw._build_support import (
    get_excluded_package_data_patterns,
    get_package_data_patterns,
    is_browser_wheel_build,
    remove_browser_wheels_from_build_dir,
    sync_browser_wheel,
)


class build_py(_build_py):
    def run(self):
        sync_browser_wheel()
        remove_browser_wheels_from_build_dir(self.build_lib)
        super().run()


class sdist(_sdist):
    def run(self):
        sync_browser_wheel()
        super().run()


INCLUDE_BROWSER_WHEEL = not is_browser_wheel_build()


setup(
    packages=find_packages(include=["gbdraw*"]),
    package_data={"gbdraw": get_package_data_patterns(include_browser_wheel=INCLUDE_BROWSER_WHEEL)},
    exclude_package_data={"gbdraw": get_excluded_package_data_patterns(include_browser_wheel=INCLUDE_BROWSER_WHEEL)},
    include_package_data=True,
    cmdclass={"build_py": build_py, "sdist": sdist},
)
