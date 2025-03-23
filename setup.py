from setuptools import setup, find_packages

setup(
    name="hap2num",
    version="0.1.0",
    author="Félicien Akohoue",
    author_email="akohoue.f@gmail.com",
    description="A package for converting haplotype genotype data to numeric format",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/FAkohoue/hap2num",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "tqdm"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)