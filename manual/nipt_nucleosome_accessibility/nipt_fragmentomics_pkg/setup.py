from setuptools import setup, find_packages

setup(
    name="nipt_fragmentomics",
    version="0.1.0",
    description="NIPT cfDNA fragmentomics pipeline: CNV, WPS, Fetal Fraction",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.24",
        "pandas>=2.0",
        "pysam>=0.21",
        "pyarrow>=12.0",
        "statsmodels>=0.14",
        "scipy>=1.11",
        "plotly>=5.18",
        "pyBigWig>=0.3",
    ],
    extras_require={
        "dev": ["pytest>=7.4", "pytest-cov"],
    },
    entry_points={
        "console_scripts": [
            "nipt-run=run:main",
        ],
    },
)
