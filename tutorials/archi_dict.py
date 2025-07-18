archi = {
    "nb_phy": [8,22], # number of phytomers on the main stem
    "nb_short_phy": [2,4],
    "short_phy_height": [2,5],
    
    # Stem
    # height=1.2*max(height_canopy), # potential plant height
    "stem_q": [0.9,1.1], # parameter for ligule height geometric distribution along axis
    "diam_base": 2.5, # stem base diameter cm
    "diam_top": 1.5, # stem top diameter cm

    # Leaf area distribution along the stem  
    # "leaf_area":[], # potential plant leaf area
    "rmax": [0.5,0.9], # relative position of largest leaf on the stem
    "skew": 0.005, # skewness for leaf area distribution along axis

    # blade area
    "wl": [0.1,0.13], # leaf blade width-to-length ratio 
    "klig": 0.6, # parameter for leaf blade shape
    "swmax": 0.55, # relative position of maximal blade width
    "f1": 0.64, # parameter for leaf blade shape
    "f2": 0.92, # parameter for leaf blade shape

    # blade curvature
    "insertion_angle": 35, # leaf blade insertion angle
    "scurv": 0.7, #  relative position of inflexion point
    "curvature": 120, # leaf blade insertion-to-tip angle
    "phyllotactic_angle": 137.5, # phyllotactic angle
    "phyllotactic_deviation": 0, # half-deviation to phyllotactic angle

    # Development
    "phyllochron": [30,60], # phyllochron, i.e. stem element appearance rate
    "plastochron": [30,60], # plastochron, i.e. leaf blade appearance rate

    # Senescence 
    # leaf_lifespan = [lifespan_early, lifespan], # leaf lifespan from appearance

    # Tillering
    "nb_tillers": 0, # number of tillers
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller
    "tiller_angle": 30, # degree
    "reduction_factor": 1, # reduction factor between tillers of consecutive order

    "plant_orientation": 0 # degree
}