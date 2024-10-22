from setuptools import setup

from treesort.version import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    install_requires=[
        'scipy>=1.7.0',
        'biopython>=1.67',
        'dendropy>=4.5.0',
        'phylo-treetime>=0.9.4',
        'matplotlib'
    ],
    name="TreeSort",
    version=__version__,
    author="Alexey Markin",
    author_email="alex.markin57@gmail.com",
    license='MIT',
    description="Virus reassortment inference software."
                "Infers both recent and ancestral reassortment and uses flexible molecular clock constraints.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/flu-crew/TreeSort",
    packages=["treesort"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={"console_scripts": ["treesort=treesort.cli:run_treesort_cli"]},
    py_modules=["treesort"],
)
