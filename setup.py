from setuptools import setup, find_packages

setup(
    name="valyte",
    version="0.1.7",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "pymatgen",
    ],
    include_package_data=True,
    package_data={
        "valyte": ["*.png"],
    },
    entry_points={
        "console_scripts": [
            "valyte=valyte.cli:main",
        ],
    },
    author="Nikhil",
    author_email="nikhil@example.com",
    description="A comprehensive CLI tool for VASP pre-processing (Supercells, K-points) and post-processing (DOS, Band Structure plotting)",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/nikyadav002/Valyte-Project",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
