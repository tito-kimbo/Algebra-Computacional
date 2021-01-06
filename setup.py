import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="python-alcp-2020-mrv", # Replace with your own username
    version="0.0.1",
    author="Eduardo Rivero, Alberto Maurel, Pablo Villalobos",
    description="A python-based computer algebra system",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tito-kimbo/Algebra-Computacional",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=[],
    python_requires='>=3.6',
)
