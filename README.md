# ArchiCrop

ArchiCrop is a 3D+t parametric plant model for cereals that gives architecture to 1D+t crop-scale data.

ArchCrop is a model where:
- plant architecture and development are defined by parameters;
- known allometrics laws enables to define a potential plant with relatively few parameters;
- plant growth is globally constrained by dynamics at crop scale, from crop model outputs for example. 

Our approach is inspired by the concept of equifinality, which suggests that an exact same outcome can be achieved starting from distinct initial conditions. 
Architectural parameters constitute degrees of freedom with respect to the crop-scale data, i.e. all the 3D plants simulated through time are equivalent regarding the crop dynamics. 
Local differential growth of each organ is calculated from a global constraint, that integrates environmental stresses. 
Based on these constraints we are able to generate a morphospace of growing architectured plants that follow exactly the crop dynamic constraints. 

The explicit spatial arrangement of crops is considered in ArchiCrop, thanks to the individual-based approach. 


## Installation
- Clone ArchiCrop repo on your own device : 
```git clone https://github.com/orianebro/ArchiCrop.git```
- Move into the repo : ```cd ./ArchiCrop/```
- Create appropriate conda environment :
  ```conda env create -f archicrop.yml ```
- Install in user mode :
  ```python setup.py install```
  or install in develop mode :  
  ```python setup.py develop```


## Tutorials
Link to notebooks


## Status 
[![Actions Status][actions-badge]][actions-link]
[![Documentation Status][rtd-badge]][rtd-link]

[![PyPI version][pypi-version]][pypi-link]
[![Conda-Forge][conda-badge]][conda-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

[![GitHub Discussion][github-discussions-badge]][github-discussions-link]

<!-- SPHINX-START -->

<!-- prettier-ignore-start -->
[actions-badge]:            https://github.com/orianebro/ArchiCrop/workflows/CI/badge.svg
[actions-link]:             https://github.com/orianebro/ArchiCrop/actions
[conda-badge]:              https://img.shields.io/conda/vn/openalea3/ArchiCrop
[conda-link]:               https://github.com/conda-forge/ArchiCrop-feedstock
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/orianebro/ArchiCrop/discussions
[pypi-link]:                https://pypi.org/project/ArchiCrop/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/ArchiCrop
[pypi-version]:             https://img.shields.io/pypi/v/ArchiCrop2
[rtd-badge]:                https://readthedocs.org/projects/ArchiCrop2/badge/?version=latest
[rtd-link]:                 https://ArchiCrop2.readthedocs.io/en/latest/?badge=latest

<!-- prettier-ignore-end -->
