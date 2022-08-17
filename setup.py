from pathlib import Path

from setuptools import find_packages, setup

with open("README.rst") as f:
    readme = f.read()

with Path("acceltools").joinpath("__init__.py").open("r") as f:
    _version = f.readline().split()[2].replace("'", "").replace('"', "")

setup(
    name="acceltools",
    version=_version,
    description="tools for ACCeL",
    long_description=readme,
    author="Keisuke Fukaya",
    author_email="kfukaya@pu-toyama.ac.jp",
    url="https://github.com/kfchem/acceltools",
    license="MIT License",
    install_requires=["accel>=0.1.0", "pandas", "scipy", "matplotlib", "python-docx"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    packages=find_packages(exclude=("tests", "docs", "application")),
    include_package_data=True,
)
