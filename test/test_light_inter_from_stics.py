import pandas as pd
from openalea.astk.sky_irradiance import sky_irradiance
from openalea.astk.sky_sources import sky_sources, caribu_light_sources
from openalea.archicrop.ltfs import illuminate
from openalea.archicrop.display import build_scene
from openalea.plantgl.all import *


fn = 'climsorj.meteo'
def meteo_day():
    names=['station', 'year', 'month', 'day', 'julian', 'min_temp', 'max_temp', 'rad', 'Penman PET', 'rainfall', 'wind', 'pressure', 'CO2']
    df = pd.read_csv(fn,  header=None, sep='\s+', names=names)
    df["daydate"] = pd.to_datetime(df[["year", "month", "day"]])
    return df

def create_plane(size=1):
    half = size / 2

    points = [Vector3(-half, -half, 0),
            Vector3(half, -half, 0),
            Vector3(half, half, 0),
            Vector3(-half, half, 0)]

    indices = [[0, 1, 2, 3]]
    faceset = FaceSet(points, indices)
    return Shape(faceset) #, Material(Color3(*color), transparency=0.5))

def test_day_to_hour():

    df = meteo_day()
    location ={
    'longitude': 3.87,
    'latitude': 45,
    'altitude': 56,
    'timezone': 'Europe/Paris'}

    scene = Scene()
    scene.add(Shape(Translated((0,0,50), Oriented(Vector3((1,0,0)), Vector3((0,1,0)),create_plane().geometry))))

    nrj = []

    for row in df.itertuples():
        irr = sky_irradiance(daydate=row.daydate, day_ghi=row.rad, **location)
        sun, sky = sky_sources(sky_type='blended', sky_irradiance=irr, scale='global')
        lights = caribu_light_sources(sun, sky)
        print(lights)
        # scene, labels = build_scene(mtg, (0,0,0), senescence=False)
        cs, raw, agg = illuminate(scene, light=lights, domain=((-1,-1),(1,1)))
        nrj.append(agg['Eabs']*agg['area'])

    print(nrj)


    
