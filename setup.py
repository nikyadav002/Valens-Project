from setuptools import setup, find_packages

setup(
    name="valens",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "pymatgen",
    ],
    entry_points={
        "console_scripts": [
            "valens=valens.cli:main",
        ],
    },
    author="Nikhil",
    description="A CLI tool for VASP post-processing (DOS plotting)",
)
