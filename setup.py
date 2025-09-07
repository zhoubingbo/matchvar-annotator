#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import os

# Read README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "MATCHVAR Annotator - Functional annotation and analysis of genomic variants"

# Read requirements file
def read_requirements():
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    if os.path.exists(requirements_path):
        with open(requirements_path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return [
        'pandas>=1.3.0',
        'psutil>=5.8.0',
        'numpy>=1.21.0'
    ]

setup(
    name="matchvar-annotator",
    version="1.0.0",
    author="Bingbo Zhou",
    author_email="zhoubingbo@hotmail.com",
    description="MATCHVAR Annotator - Functional annotation and analysis of genomic variants",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/zhoubingbo/matchvar-annotator",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.0',
            'black>=21.0',
            'flake8>=3.8',
            'mypy>=0.800',
        ],
        'docs': [
            'sphinx>=4.0',
            'sphinx-rtd-theme>=0.5',
        ]
    },
    entry_points={
        'console_scripts': [
            'matchvar-annotator=matchvar_annotator.cli:main',
            'matchvar-table=matchvar_annotator.table_matchvar:main',
            'matchvar-convert=matchvar_annotator.convert2matchvar:main',
            'matchvar-coding=matchvar_annotator.coding_change:main',
            'matchvar-index=matchvar_annotator.build_tabix_indexes:main',
            'matchvar-db=matchvar_annotator.db_cli:main',
        ],
    },
    include_package_data=True,
    package_data={
        'matchvar_annotator': [
            'data/*',
            'templates/*',
            '*.py',
        ],
    },
    zip_safe=False,
    keywords="bioinformatics, variant annotation, genomics, MATCHVAR",
    project_urls={
        'Bug Reports': 'https://github.com/zhoubingbo/matchvar-annotator/issues',
        'Source': 'https://github.com/zhoubingbo/matchvar-annotator',
        'Documentation': 'https://matchvar-annotator.readthedocs.io/',
    },
)
