# ipython --gui=qt
from __future__ import annotations

import numpy as np
import pandas as pd

from . import sky_sources as skys
from .cereals import build_shoot, leaf_azimuth, parametric_leaf
from .display import build_scene, display_scene
from .geometry import form_factor
from .simple_maize import bell_shaped_dist, geometric_dist
from .stand import agronomic_plot

# import light_it as ltfs


# generation of a 3D plant from descritive parameters
stem_radius = 0.5
insertion_heights = (10, 20)
leaf_lengths = (10, 10)
leaf_areas = (10, 10)
# type ?parametric_leaf for parameter siginification
a_leaf = parametric_leaf(
    nb_segment=10, insertion_angle=50, scurv=0.5, curvature=50, alpha=-2.3
)
leaf_shapes = [a_leaf for l in leaf_lengths]
# type ?leaf_azimuths for parameter siginification
leaf_azimuths = leaf_azimuth(
    size=len(leaf_lengths),
    phyllotactic_angle=180,
    phyllotactic_deviation=15,
    plant_orientation=0,
    spiral=False,
)
#
shoot, g = build_shoot(
    stem_radius=stem_radius,
    insertion_heights=insertion_heights,
    leaf_lengths=leaf_lengths,
    leaf_areas=leaf_areas,
    leaf_shapes=leaf_shapes,
    leaf_azimuths=leaf_azimuths,
)
scene, nump = build_scene(g)
display_scene(scene)

# some plausible value for a maize plant
stem_radius = 1
leaf_areas = bell_shaped_dist(plant_area=10000)
a_leaf = parametric_leaf(
    nb_segment=10, insertion_angle=50, scurv=0.5, curvature=50, alpha=-2.3
)
ff = form_factor(a_leaf)
leaf_lengths = np.sqrt(np.array(leaf_areas) / 0.1 / ff)
insertion_heights = np.cumsum(geometric_dist(height=200))
leaf_shapes = [a_leaf for l in leaf_lengths]
leaf_azimuths = leaf_azimuth(
    size=len(leaf_lengths),
    phyllotactic_angle=180,
    phyllotactic_deviation=15,
    plant_orientation=0,
    spiral=False,
)
shoot, g = build_shoot(
    stem_radius=stem_radius,
    insertion_heights=insertion_heights,
    leaf_lengths=leaf_lengths,
    leaf_areas=leaf_areas,
    leaf_shapes=leaf_shapes,
    leaf_azimuths=leaf_azimuths,
)
scene, nump = build_scene(g)
display_scene(scene)


# some realistic values for a wheat plant
df = pd.read_csv("wheat.csv")
stem_radius = 0.25
a_leaf = parametric_leaf(
    nb_segment=10, insertion_angle=30, scurv=0.5, curvature=20, alpha=-1.5
)
ff = form_factor(a_leaf)
leaf_areas = df.L_blade * df.W_blade * ff
leaf_lengths = df.L_blade
insertion_heights = df.H_blade
leaf_shapes = [a_leaf for l in leaf_lengths]
leaf_azimuths = leaf_azimuth(
    size=len(leaf_lengths),
    phyllotactic_angle=600,
    phyllotactic_deviation=10,
    plant_orientation=0,
    spiral=True,
)


# generate x,y position for a stand
nplants, positions, domain, domain_area, unit = agronomic_plot(1, 1, 10, 10, 0.75)
shoot, g = build_shoot(
    stem_radius=stem_radius,
    insertion_heights=insertion_heights,
    leaf_lengths=leaf_lengths,
    leaf_areas=leaf_areas,
    leaf_shapes=leaf_shapes,
    leaf_azimuths=leaf_azimuths,
)
plants = [g for i in range(nplants)]
scene, nump = build_scene(plants, positions)
display_scene(scene)

# vertical light interception
# cs, ei, df = ltfs.illuminate(scene, scene_unit='cm')
# cs.plot(ei)


# diffuse light interception
sources = skys.sky_sources()
# cs, ei, df = ltfs.illuminate(scene, light=sources, scene_unit='cm')
# cs.plot(ei)

# get score per plant
# def score(res):
#     return pandas.Series({'ei':(res.Ei*res.area).sum() / res.area.sum(),
#                               'area': res.area.sum()})
# df['nump']=nump
# df.groupby('nump').apply(score)
