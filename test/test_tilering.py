from openalea.archicrop.generator_bis import cereals
from openalea.archicrop.cereals import build_shoot

def test1():

    
    nb_phy=12
    height=200
    leaf_area=4000
    wl=0.12
    diam_base=2
    diam_top=1
    insertion_angle=30
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
    phyllotactic_deviation=20

    xxx= build_shoot(
    nb_phy,
    height,
    leaf_area,
    wl,
    diam_base,
    diam_top,
    insertion_angle,
    scurv,
    curvature,
    klig, swmax, f1, f2,
    stem_q,
    rmax,
    skew,
    phyllotactic_angle,
    phyllotactic_deviation)

    return xxx