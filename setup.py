#!/usr/bin/python3

from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(name='provis',
      version='0.0.3',
      description='Protein Visualization Library in Python',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/czirjakkethz/provis',
      author='Kristof Czirjak',
      author_email='czirjakk@student.ethz.ch',
      license='MIT',
      packages=find_packages(),
      install_requires=['biopython', 'trimesh', 'pyvista', 'biopandas', 'torch', 'pyvtk', 'open3d'],
      keywords=['python', 'protein', 'visualization', 'pdb'],
      zip_safe=False,

)
