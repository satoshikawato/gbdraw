from pathlib import Path

from setuptools import find_packages, setup
from setuptools.command.build_py import build_py as _build_py


class build_py(_build_py):
    def run(self):
        super().run()

        # Keep build/lib in sync with the source tree so removed browser wheels
        # do not linger in subsequent wheel builds.
        source_web_root = Path("gbdraw") / "web"
        build_web_root = Path(self.build_lib) / "gbdraw" / "web"
        source_wheels = {path.name for path in source_web_root.glob("*.whl")}
        for wheel_path in build_web_root.glob("*.whl"):
            if wheel_path.name not in source_wheels:
                wheel_path.unlink()

setup(
    name='gbdraw',
    version='0.9.1',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'svgwrite',
        'biopython',
        'fonttools',
        'bcbio-gff'
    ],
    extras_require={
        'dev': [
            'cairosvg',
            'Pillow',
            'pytest>=7.0',
            'pytest-cov>=4.0',
            'pytest-timeout>=2.0',
            'setuptools>=61.0',
            'wheel',
        ],
        'export': ['cairosvg'],
    },
    package_data={
        'gbdraw': [
            'data/color_palettes.toml', 
            'data/config.toml', 
            'data/*.ttf',
            'web/*.whl',
            'web/index.html',
            'web/open-source-notices.html',
            'web/js/*.js',
            'web/js/app/*.js',
            'web/js/app/*/*.js',
            'web/js/workers/*.js',
            'web/js/services/*.js',
            'web/js/utils/*.js',
            'web/presets/*.txt',
            'web/presets/*.tsv',
            'web/vendor/fonts/*.ttf',
            'web/vendor/*/*.css',
            'web/vendor/*/*.js',
            'web/vendor/*/*.json',
            'web/vendor/*/*.svg',
            'web/vendor/*/*.ttf',
            'web/vendor/*/*.whl',
            'web/vendor/*/*.woff',
            'web/vendor/*/*.woff2',
            'web/vendor/*/*/*.css',
            'web/vendor/*/*/*.js',
            'web/vendor/*/*/*.json',
            'web/vendor/*/*/*.mjs',
            'web/vendor/*/*/*.svg',
            'web/vendor/*/*/*.ttf',
            'web/vendor/*/*/*.whl',
            'web/vendor/*/*/*.woff',
            'web/vendor/*/*/*.woff2',
            'web/vendor/*/*/*.zip',
            'web/vendor/*/*/*.wasm',
            'web/vendor/*/*/*/*.js',
            'web/vendor/*/*/*/*.json',
            'web/vendor/*/*/*/*.mjs',
            'web/vendor/*/*/*/*.whl',
            'web/vendor/*/*/*/*.zip',
            'web/vendor/*/*/*/*.wasm',
            'web/wasm/*/*.md',
            'web/wasm/*/*.wasm',
        ]
    },
    include_package_data=True,
    python_requires='>=3.10',
    author='Satoshi Kawato',
    author_email='kawato@kaiyodai.ac.jp',
    description='a genome diagram generator for microbes and organelles',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/satoshikawato/gbdraw/',
    project_urls={
        "Web App": "https://gbdraw.app/",
        "Source": "https://github.com/satoshikawato/gbdraw/",
    },
    entry_points={'console_scripts': ['gbdraw = gbdraw.cli:main',] },
    cmdclass={"build_py": build_py},
)
