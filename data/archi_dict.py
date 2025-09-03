from __future__ import annotations

archi_maize = {
    # "nb_phy": [8,25], # number of phytomers on the main stem
    "nb_short_phy": [2,5],
    "short_phy_height": [2,5],
    
    # Stem
    "stem_q": [0.8,1.2], # parameter for ligule height geometric distribution along axis
    "diam_base": 2.2, # stem base diameter cm
    "diam_top": 1.2, # stem top diameter cm

    # Leaf area distribution along the stem  
    "rmax": [0.5,0.9], # relative position of largest leaf on the stem
    "skew": [0.0001,0.01], # skewness for leaf area distribution along axis 

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
    "phyllochron": [20,55], # phyllochron, i.e. stem element appearance rate
    "plastochron": [30,60], # plastochron, i.e. leaf blade appearance rate
    "stem_duration": 1.6,
    "leaf_duration": 1.6,

    # Tillering
    "nb_tillers": 0, # number of tillers
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller
    "tiller_angle": 30, # degree
    "reduction_factor": 1, # reduction factor between tillers of consecutive order

    "plant_orientation": 0 # degree
}


archi_sorghum = {
    # "nb_phy": [8,25], # number of phytomers on the main stem : [10,40] (Ndiaye et al., 2021; Lafarge and Tardieu, 2002; Clerget, 2008; Ganeme et al., 2022)
    "nb_short_phy": [2,5],
    "short_phy_height": [2,5],

    # Stem
    "stem_q": [0.9,1.1], # parameter for ligule height distribution along axis : [1.1] (Kaitaniemi et al., 1999) 
    "diam_base": 2.2, # stem base diameter : [2.2] (Ndiaye et al., 2021)
    "diam_top": 1.2, # stem top diameter: [1.2] (Ndiaye et al., 2021)

    # Leaf area distribution along the stem
    "rmax": [0.5,0.9], # parameter for leaf area distribution along axis : [0.6,0.8] (Kaitaniemi et al., 1999; Welcker et al., )
    "skew": [0.0001,0.001], # parameter for leaf area distribution along axis : [0.0005,0.1] (Kaitaniemi et al., 1999; Welcker et al., )
    
    # blade area 
    "wl": [0.1,0.12], # leaf blade width-to-length ratio : [0.1, 0.12] ()
    "klig": 0.6, # parameter for leaf blade shape
    "swmax": 0.55, # parameter for leaf blade shape
    "f1": 0.64, # parameter for leaf blade shape
    "f2": 0.92, # parameter for leaf blade shape

    # Leaf blade position in space
    "insertion_angle": 30, # leaf blade insertion angle : [10,50] (Truong et al., 2015; Kaitaniemi et al., 1999)
    "scurv": 0.7, # leaf blade relative inflexion point : [0.6, 0.8] ()
    "curvature": 90, # leaf blade insertion-to-tip angle : [45, 135] (Kaitaniemi et al., 1999)
    "phyllotactic_angle": 180, # phyllotactic angle : [180] (Davis et al., 2024)
    "phyllotactic_deviation": 10, # half-deviation to phyllotactic angle : [0,90] (Davis et al., 2024)

    # Development
    "phyllochron": [20,70], # phyllochron, i.e. stem element appearance rate : [40,65 then x1.6-2.5] (Clerget, 2008)
    # "plastochron": [40,65], # plastochron, i.e. leaf blade appearance rate : [34,46 then 80-93] (Rami Kumar et al., 2009)
    "stem_duration": [1.6,1.8],
    "leaf_duration": [1.6,1.8],

    # Tillering
    "nb_tillers": 0, # number of tillers : [0,6] (Lafarge et al., 2002)
    "tiller_angle": 20,
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller : [] ()
    "reduction_factor": 1 # reduction factor between tillers of consecutive order : [0.8,1] ()
}


archi_sorghum_angles = {
    # "nb_phy": 12, # number of phytomers on the main stem : [10,40] (Ndiaye et al., 2021; Lafarge and Tardieu, 2002; Clerget, 2008; Ganeme et al., 2022)
    "nb_short_phy": 4,
    "short_phy_height": 3,

    # Stem
    "stem_q": 0.9, # parameter for ligule height distribution along axis : [1.1] (Kaitaniemi et al., 1999) 
    "diam_base": 2.2, # stem base diameter : [2.2] (Ndiaye et al., 2021)
    "diam_top": 1.2, # stem top diameter: [1.2] (Ndiaye et al., 2021)

    # Leaf area distribution along the stem
    "rmax": 0.67, # parameter for leaf area distribution along axis : [0.6,0.8] (Kaitaniemi et al., 1999; Welcker et al., )
    "skew": 0.0001, # parameter for leaf area distribution along axis : [0.0005,0.1] (Kaitaniemi et al., 1999; Welcker et al., )
    
    # blade area 
    "wl": 0.12, # leaf blade width-to-length ratio : [0.1, 0.12] ()
    "klig": 0.6, # parameter for leaf blade shape
    "swmax": 0.55, # parameter for leaf blade shape
    "f1": 0.64, # parameter for leaf blade shape
    "f2": 0.92, # parameter for leaf blade shape

    # Leaf blade position in space
    "insertion_angle": [10,50], # leaf blade insertion angle : [10,50] (Truong et al., 2015; Kaitaniemi et al., 1999)
    "scurv": [0.6,0.8], # leaf blade relative inflexion point : [0.6, 0.8] ()
    "curvature": [45,135], # leaf blade insertion-to-tip angle : [45, 135] (Kaitaniemi et al., 1999)
    "phyllotactic_angle": [130,180], # phyllotactic angle : [135.7;180] (Davis et al., 2024)
    "phyllotactic_deviation": 0, # half-deviation to phyllotactic angle : [0,90] (Davis et al., 2024)

    # Development
    "phyllochron": 54, # phyllochron, i.e. stem element appearance rate : [40,65 then x1.6-2.5] (Clerget, 2008)
    "plastochron": 54, # plastochron, i.e. leaf blade appearance rate : [34,46 then 80-93] (Rami Kumar et al., 2009)
    "stem_duration": 1.6,
    "leaf_duration": 1.6,

    # Tillering
    "nb_tillers": 0, # number of tillers : [0,6] (Lafarge et al., 2002)
    "tiller_angle": 20,
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller : [] ()
    "reduction_factor": 1 # reduction factor between tillers of consecutive order : [] ()
}