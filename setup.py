#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re

from setuptools import find_packages, setup

HERE = os.path.abspath(os.path.dirname(__file__))


def read_file(*parts):
    """Read a UTF-8 text file from the project root."""
    path = os.path.join(HERE, *parts)
    if not os.path.exists(path):
        return ""
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def get_version():
    """Read the version from matchvar_annotator/__init__.py (single source of truth)."""
    init_py = read_file("matchvar_annotator", "__init__.py")
    match = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', init_py, re.MULTILINE)
    if not match:
        raise RuntimeError("Cannot find __version__ in matchvar_annotator/__init__.py")
    return match.group(1)


def get_requirements():
    """Parse requirements.txt, ignoring comments and blank lines."""
    requirements = []
    for line in read_file("requirements.txt").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            requirements.append(line)
    return requirements or [
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "psutil>=5.8.0",
    ]


DESCRIPTION = "MATCHVAR Annotator - Functional annotation and analysis of genomic variants"

setup(
    name="matchvar-annotator",
    version=get_version(),
    author="Bingbo Zhou",
    author_email="zhoubingbo@hotmail.com",
    description=DESCRIPTION,
    long_description=read_file("README.md") or DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/zhoubingbo/matchvar-annotator",
    license="MIT",
    license_files=("LICENSE",),
    packages=find_packages(exclude=("tests", "tests.*")),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Environment :: Console",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=get_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.800",
            "build>=1.0.0",
            "twine>=4.0.0",
        ],
        "docs": [
            "sphinx>=4.0",
            "sphinx-rtd-theme>=0.5",
        ],
        "bb": [
            "pyBigWig>=0.3.18",
        ],
    },
    entry_points={
        "console_scripts": [
            "matchvar-annotator=matchvar_annotator.cli:main",
            "matchvar-table=matchvar_annotator.table_matchvar:main",
            "matchvar-convert=matchvar_annotator.convert2matchvar:main",
            "matchvar-coding=matchvar_annotator.coding_change:main",
            "matchvar-index=matchvar_annotator.build_tabix_indexes:main",
            "matchvar-db=matchvar_annotator.db_cli:main",
            "matchvar-clinvardb=matchvar_annotator.vcf2clinvar_db:main",
            "matchvar-dbnsfpdb=matchvar_annotator.dbnsfp2matchvar_db:main",
            "matchvar-pipeline=matchvar_annotator.simulate_annotate_pipeline:main",
            "matchvar-ghgvs=matchvar_annotator.ghgvs:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords="bioinformatics, variant annotation, genomics, MATCHVAR",
    project_urls={
        "Bug Reports": "https://github.com/zhoubingbo/matchvar-annotator/issues",
        "Source": "https://github.com/zhoubingbo/matchvar-annotator",
        "Documentation": "https://matchvar-annotator.readthedocs.io/",
        "Changelog": "https://github.com/zhoubingbo/matchvar-annotator/blob/main/CHANGELOG.md",
    },
)
