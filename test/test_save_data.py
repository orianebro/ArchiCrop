from __future__ import annotations

import pandas as pd
import xarray as xr

dates = xr.date_range("2000-01-01", periods=5)

params = {
    0: {"a" : 1, "b" : 2, "c" : 3, "d" : 4, "e" : 5},
    1: {"a" : 6, "b" : 7, "c" : 8, "d" : 9, "e" : 10},
    2: {"a" : 11, "b" : 12, "c" : 13, "d" : 14, "e" : 15}
}

daily_dyn = {
    "tps_ther": [20, 21, 22, 23, 24],
    "lai": [100, 101, 102, 103, 104]
}

pot_la = {
    0: [1, 2, 3, 4, 5],
    1: [2, 3, 4, 5, 6],
    2: [3, 4, 5, 6, 7]
}

pot_h = {
    0: [10, 11, 12, 13, 14],
    1: [11, 12, 13, 14, 15],
    2: [12, 13, 14, 15, 16]
}

ds = xr.Dataset(
    data_vars=dict(
        STICS_tps_ther = (["time"], daily_dyn["tps_ther"]),
        STICS_lai = (["time"], daily_dyn["lai"]),
        pot_lat = (["id", "time"], pd.DataFrame.from_dict(pot_la, orient='index', columns=dates)),
        pot_h = (["id", "time"], pd.DataFrame.from_dict(pot_h, orient='index', columns=dates)),
    ),
    coords=dict(
        id = [0,1,2],
        time = dates
    )
)
