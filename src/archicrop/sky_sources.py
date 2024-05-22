from __future__ import annotations

import numpy as np
import openalea.plantgl.all as pgl
import pandas as pd
from alinea.astk.meteorology.sky_irradiance import sun_position
from alinea.astk.sun_and_sky import sun_sources as ss

from .colormap import jet_colors

# from alinea.astk.sun_and_sky import sky_sources, sun_sky_sources


def cartesian(elevation, azimuth, radius=1):
    """Cartesian coordinates from elevation and azimuth angles"""
    zenith = np.radians(90 - elevation)
    azimuth = np.radians(azimuth)
    x = radius * np.sin(zenith) * np.cos(azimuth)
    y = radius * np.sin(zenith) * np.sin(azimuth)
    z = radius * np.cos(zenith)
    return zip(x, y, z)


def show_sources(sun, north=0, scene=None, distance=5, radius=0.1):
    sc = pgl.Scene()
    if scene is not None:
        sc.add(scene)
    elevation, azimuth, luminance = (np.array(x) for x in sun)
    az = -(azimuth + north)
    colors = jet_colors(luminance)
    pgl_colors = [pgl.Material(pgl.Color3(r, g, b)) for r, g, b in colors]
    sph = pgl.Sphere(radius)
    pos = cartesian(elevation, az, distance)
    pgl_pos = [pgl.Vector3(*p) for p in pos]
    for i, vect in enumerate(pgl_pos):
        shape = pgl.Translated(vect, sph)
        sc.add(pgl.Shape(shape, pgl_colors[i]))
    pgl.Viewer.display(sc)
    return sc


def daylength(daydate, **location):
    sun = sun_position(daydate=daydate, **location)
    return len(sun)


def add_sources(sun, sky):
    el, az, lum = (list(x) for x in sky)
    elevation, azimuth, luminance = (list(x) for x in sun)
    return el + elevation, az + azimuth, lum + luminance


def sun_sources(irradiance=1, *args, **kwds):
    irrad = irradiance
    if irradiance == "clear_sky":
        irrad = None
    return ss(irrad, *args, **kwds)


def read_meteo_mpt(when="winter"):
    if when == "winter":
        path = "incoming_radiation_ZA13.csv"
    else:
        path = "incoming_radiation_ZB13.csv"
    df = pd.read_csv(
        path,
        index_col=0,
        parse_dates=True,
        date_parser=lambda col: pd.to_datetime(col, utc=True),
    )
    df.index = df.index.tz_convert("Europe/Paris")

    return df


def montpellier_spring_2013():
    return read_meteo_mpt("spring")


def montpellier_winter_2013():
    return read_meteo_mpt("winter")


def actual_irradiance(daydate, db):
    when = pd.date_range(start=daydate, freq="H", periods=25, tz="Europe/Paris")
    return db.loc[when, "ghi"]
