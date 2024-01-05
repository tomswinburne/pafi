import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PAFI",
    version="0.9.0",
    author="TD Swinburne",
    author_email="thomas.swinburne@cnrs.fr",
    description="PAFI: MD evaluation of free energy barriers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tomswinburne/pafi",
    packages=setuptools.find_packages(),
    install_requires=['mpi4py>=3.1.4','scipy>=1.10','numpy>=1.17','pandas>=1.5.0','plotext>=5.2.8','tqdm>=4.50.0'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
)