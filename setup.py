from setuptools import setup

with open("README.md","r") as fh: 
    long_description = fh.read()

setup(
    name="nais-processor",
    version="0.0.12",
    description='Code to process ion spectrometer data files',
    py_modules=["nais_processor"],
    package_dir={'':'src'},
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires = [
        "pandas >= 1.1.0",
        "numpy >= 1.19.0",
        "matplotlib >= 3.3.4",
        "pyyaml >= 5.4.1",
        "tinydb >= 4.7.0",
        "aerosol-functions >= 0.0.2"
    ],
    url="https://github.com/jlpl/nais-processor",
    author="Janne Lampilahti",
    author_email="janne.lampilahti@helsinki.fi",
)
