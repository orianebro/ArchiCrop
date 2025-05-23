from __future__ import annotations

import numpy as np

from .cereals_leaf import parametric_leaf
from .generator import cereals as cereals_generator
from .plant_design import blade_dimension, leaf_azimuth, stem_dimension
from .plant_shape import bell_shaped_dist, geometric_dist


def build_shoot(
    nb_phy,
    height,
    leaf_area,
    wl,
    diam_base,
    diam_top,
    insertion_angle,
    scurv,
    curvature,
    # alpha,
    klig, swmax, f1, f2,
    stem_q,
    rmax,
    skew,
    phyllotactic_angle,
    phyllotactic_deviation,
    plant_orientation=45,
    spiral=True,
    nb_short_phy=4
):
    """create a shoot

    Args:
        stem_radius: (float) the stem radius
        insertion_heights: list of each leaf insertion height
        leaf_lengths: list of each leaf length (blade length)
        leaf_areas: list of each blade area
        collar_visible: list of each collar height or True if the collar is visible and False if it is not
        leaf_shapes: list of each leaf shape, if it is not known write None
        leaf_azimuths: list of each leaf azimuth, if it is not known write None

    Returns:
        shoot:

    """

    ranks = range(1, nb_phy + 1)
    ntop = max(ranks) - np.array(ranks) + 1

    ## Stem
    # internode lengths
    #nb_young_phy = int(
    #    round((nb_phy - 1.95) / 1.84 / 1.3)
    #)  # Lejeune and Bernier formula + col

    # nb_young_phy = int(
    #     round((nb_phy - 1.95) / 1.84 / 1.3)
    # )
    # nb_short_phy = 4
    short_phy_height = 2

    pseudostem_height = nb_short_phy * short_phy_height

    pseudostem = np.array([short_phy_height*i for i in range(1, nb_short_phy+1)])
    stem = np.array(geometric_dist(height, nb_phy-nb_short_phy, q=stem_q, u0=pseudostem_height))
    insertion_heights = np.concatenate((pseudostem, stem), axis=0)
    # stem = np.array([pseudostem_height+i*(height - pseudostem_height)/(nb_phy - nb_young_phy) for i in range(nb_young_phy+1,nb_phy+1)])
    # insertion_heights = np.array(geometric_dist(height, nb_phy, q=stem_q)) #, u0=young_phy_height))
    # insertion_heights = np.array(collar_heights_kaitaniemi(height, nb_phy))

    # stem diameters
    # stem_diameters = [diam_base] * nb_young_phy + np.linspace(
    #     diam_base, diam_top, nb_phy - nb_young_phy
    # ).tolist()
    stem_diameters = np.linspace(diam_base, diam_top, nb_phy).tolist()

    stem = stem_dimension(
        h_ins=insertion_heights, d_internode=stem_diameters, ntop=ntop
    )

    ## Leaves

    # leaf length repartition along axis
    # leaf_areas_stem = np.array(bell_shaped_dist(leaf_area, nb_phy, rmax, skew))
    # leaf_areas_pseudostem = np.array(bell_shaped_dist(leaf_area, nb_phy, rmax, skew))
    # leaf_areas = np.concatenate((leaf_areas_pseudostem, leaf_areas_stem), axis=0)
    leaf_areas = np.array(bell_shaped_dist(leaf_area, nb_phy, rmax, skew))


    # leaf shapes
    a_leaf = parametric_leaf(
        nb_segment=10,
        insertion_angle=insertion_angle,
        scurv=scurv,
        curvature=curvature,
        klig=klig, swmax=swmax, f1=f1, f2=f2
    )

    leaf_shapes = [a_leaf] * nb_phy 

    blades = blade_dimension(area=leaf_areas, ntop=ntop, wl=wl)

    # tck = shape_to_surface(a_leaf, wl)

    # ff = [get_form_factor(leaf) for leaf in leaf_shapes]

    # phyllotaxy
    leaf_azimuths = leaf_azimuth(
        size=nb_phy,
        phyllotactic_angle=phyllotactic_angle,
        phyllotactic_deviation=phyllotactic_deviation,
        plant_orientation=plant_orientation,
        spiral=spiral,
    )
    leaf_azimuths[1:] = np.diff(leaf_azimuths)

    ## df
    axis = blades + stem
    axis["leaf_azimuth"] = leaf_azimuths
    axis["leaf_rank"] = ranks
    axis["leaf_shape"] = [leaf_shapes[n - 1] for n in ranks]
    # df["wl"] = [wl for _ in df.leaf_rank]
    # df["leaf_tck"] = [(tck) for _ in df.leaf_rank]
    return axis, cereals_generator(plant=axis)


def shoot_at_stage(shoot, stage):
    df = shoot.loc[shoot["leaf_rank"] <= stage, :]  # noqa: PD901
    return df, cereals_generator(plant=df)
