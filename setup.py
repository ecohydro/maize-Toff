import setuptools

with open("README.md", "r") as fh:
    README_description = fh.read()

setuptools.setup(
    name="Maize_Toff",
    version="0.1",
    author="Natasha Krell",
    author_email="nkrell@ucsb.edu>",
    description="Software to model maize yield variability & tradeoffs between yield and crop failure.",
    long_description=README_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ecohydro/maize-Toff",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
