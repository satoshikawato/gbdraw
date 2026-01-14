from setuptools import setup, find_packages

setup(
    name='gbdraw',
    version='0.8.0',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'svgwrite',
        'biopython',
        'fonttools',
        'bcbio-gff'
    ],
    extras_require={
        'export': ['cairosvg'],
    },
    package_data={
        'gbdraw': [
            'data/color_palettes.toml', 
            'data/config.toml', 
            'data/*.ttf',
            'web/index.html',
            'web/js/*.js',
            'web/js/services/*.js',
            'web/js/utils/*.js',
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
)
