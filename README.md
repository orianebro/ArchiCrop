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


## Hands-on
```
from openalea.archicrop.archicrop import ArchiCrop

# Management parameters
sowing_density = 10
inter_row = 0.4

# Growth and senescence dynamics
daily_dynamics = dict of dicts, for each time step (day), a dict of values :
        - "Thermal time" (float): thermal time (in °C.day).
        - "Plant leaf area" (float): plant leaf area at a given thermal time (in cm²).
        - "Leaf area increment" (float): leaf area increment at a given thermal time (in cm²).
        - "Senescent leaf area" (float): senescent plant leaf area at a given thermal time (in cm²).
        - "Senescent leaf area increment" (float): senescent leaf area increment at a given thermal time (in cm²).
        - "Plant height" (float): plant height at a given thermal time (in cm).
        - "Height increment" (float): height increment at a given thermal time (in cm).
        - "Absorbed PAR" (float): absorbed PAR at a given thermal time (in MJ/m²)

thermal_time : list
leaf_area : list
sen_leaf_area :list
height : list

daily_dynamics = {
        i+1: {"Thermal time": thermal_time[i],
            "Phenology": 'germination' if i+1 <= emergence else 'juvenile' if emergence < i+1 <= end_juv else 'exponential' if end_juv < i+1 <= max_lai else 'repro',
            "Plant leaf area": leaf_area[i], 
            "Leaf area increment": leaf_area_incr[i], 
            "Senescent leaf area": sen_leaf_area[i],
            "Senescent leaf area increment": sen_leaf_area_incr[i],
            "Plant height": height[i], 
            "Height increment": height_incr[i], 
            "Absorbed PAR": par_abs[i]}
        for i in range(len(thermal_time))
}

# Senescence parameters
lifespan = 240
lifespan_early = 1 * lifespan

# ArchiCrop parameters
height=2*max(height)
Smax=2*max(leaf_area)
nb_phy=20
wl=0.12
diam_base=2.5 
diam_top=1.5
insertion_angle=35
scurv=0.7
curvature=100
klig=0.6
swmax=0.55
f1=0.64 
f2=0.92
stem_q=1.1
rmax=0.8
skew=0.0005
phyllotactic_angle=180
phyllotactic_deviation=10
phyllochron=30
plastochron=20
leaf_lifespan=[lifespan_early, lifespan]

# Instanciate ArchiCrop object
plant = ArchiCrop(height=height, 
                    nb_phy=nb_phy,
                    Smax=Smax,
                    wl=wl, diam_base=diam_base, diam_top=diam_top, 
                    insertion_angle=insertion_angle, scurv=scurv, curvature=curvature, 
                    klig=klig, swmax=swmax, f1=f1, f2=f2, 
                    stem_q=stem_q, rmax=rmax, skew=skew,
                    phyllotactic_angle=phyllotactic_angle,
                    phyllotactic_deviation=phyllotactic_deviation,
                    phyllochron=phyllochron, 
                    plastochron=plastochron, 
                    leaf_lifespan=leaf_lifespan,
                    increments=daily_dynamics)

# Generate a potential plant
plant.generate_potential_plant()
# plant.define_development()

# Simulate growth and senescence of this plant according to the given dynamics
growing_plant = plant.grow_plant()
```


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
