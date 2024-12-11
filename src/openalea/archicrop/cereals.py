from __future__ import annotations

import numpy as np

from .cereals_leaf import parametric_leaf
from .generator import cereals as cereals_generator
from .plant_design import blade_dimension, leaf_azimuth, stem_dimension
from .plant_shape import bell_shaped_dist, geometric_dist


def build_shoot(
    nb_phy,
    height,
    max_leaf_length,
    wl,
    diam_base,
    diam_top,
    insertion_angle,
    scurv,
    curvature,
    alpha,
    stem_q,
    rmax,
    skew,
    phyllotactic_angle,
    phyllotactic_deviation,
    plant_orientation=45,
    spiral=True,
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
    nb_young_phy = 5
    young_phy_height = 2

    pseudostem_height = nb_young_phy * young_phy_height

    pseudostem = np.array([young_phy_height*i for i in range(1, nb_young_phy+1)])
    stem = np.array(
        geometric_dist(height - pseudostem_height, nb_phy - nb_young_phy, q=stem_q, u0=young_phy_height)
    )
    # stem = np.array([pseudostem_height+i*(height - pseudostem_height)/(nb_phy - nb_young_phy) for i in range(nb_young_phy+1,nb_phy+1)])
    insertion_heights = np.concatenate((pseudostem, stem), axis=0)
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
    leaf_lengths_stem = np.array(
        bell_shaped_dist(
            max_leaf_length=max_leaf_length, nb_phy=nb_phy - nb_young_phy, rmax=rmax, skew=skew
        )
    )
    leaf_lengths_pseudostem = np.array(
        bell_shaped_dist(
            max_leaf_length=leaf_lengths_stem[0], nb_phy=nb_young_phy, rmax=1, skew=1
        )
    )
    leaf_lengths = np.concatenate((leaf_lengths_pseudostem, leaf_lengths_stem), axis=0)
    # leaf_lengths = np.array(bell_shaped_dist(max_leaf_length=max_leaf_length, nb_phy=nb_phy, rmax=rmax, skew=skew))

    blades = blade_dimension(area=None, length=leaf_lengths, ntop=ntop, wl=wl)

    # leaf shapes
    a_leaf = parametric_leaf(
        nb_segment=10,
        insertion_angle=insertion_angle,
        scurv=scurv,
        curvature=curvature,
        alpha=alpha,
    )
    leaf_shapes = [a_leaf] * nb_phy  # or nb_phy - nb_young_phy ???

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
    df = blades.merge(stem)
    df["leaf_azimuth"] = leaf_azimuths
    df["leaf_rank"] = ranks
    df["leaf_shape"] = [leaf_shapes[n - 1] for n in df.leaf_rank]
    return df, cereals_generator(plant=df)


def shoot_at_stage(shoot, stage):
    df = shoot.loc[shoot["leaf_rank"] <= stage, :]
    return df, cereals_generator(plant=df)
