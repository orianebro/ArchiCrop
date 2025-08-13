from __future__ import annotations

import xml.etree.ElementTree as ET

import pandas as pd

from .sky_sources import meteo_day


def read_xml_file(file_xml, params):
    """
    Parses an XML file and retrieves the values of the specified parameters.

    :param file_xml: Path to the XML file.
    :param params: List of parameter names to extract.
    :return: Dictionary with parameter names as keys and extracted values.
    """
    tree = ET.parse(file_xml)
    root = tree.getroot()
    
    result = {}
    
    # Search for all 'param' and 'colonne' elements in the XML
    for elem in root.findall(".//param") + root.findall(".//colonne"):
        param_name = elem.get("nom")  # Get the name attribute
        if param_name in params:
            result[param_name] = float(elem.text.strip()) if elem.text else None
    
    return result


def read_sti_file(file_sti, density):
    """Reads a STICS mod_s*.sti output file and builds a dictionary.
    
    :param file: str, input file of STICS outputs :
        - tempeff(n) : daily efficient thermal time (°C.day)
        - laimax : canopy max LAI (m2/m2)
        - laisen(n) : senescent LAI (m2/m2)
        - hauteur : canopy height (m)
        - raint : PAR intercepted (actually, PAR absorbed) by canopy (MJ/m2)
        - trg(n) : global radiation (MJ/m2)
    :return: dict of dicts, for each time step, a dict of values from STICS outputs, converted to be used in ArchiCrop :
        - "Thermal time" (float): thermal time (in °C.day).
        - "Plant leaf area" (float): plant leaf area at a given thermal time (in cm²).
        - "Leaf area increment" (float): leaf area increment at a given thermal time (in cm²).
        - "Plant height" (float): plant height at a given thermal time (in cm).
        - "Height increment" (float): height increment at a given thermal time (in cm).
        - "Absorbed PAR" (float): absorbed PAR at a given thermal time (in MJ/m²)"""
    
    data_dict = {}
    non_zero_height_encountered = False

    with open(file_sti) as file:  # noqa: PTH123
        # Read the header line to get column names
        header = file.readline().strip().split(";")
        # Strip whitespace from column names
        stripped_header = [col.strip() for col in header if col != 'pla']

        # Find indices for date columns
        ian_idx = header.index("ian")
        mo_idx = header.index("mo")
        jo_idx = header.index("jo")

        # Initialize empty lists for each selected column in the dictionary
        data_dict = {col.strip(): [] for col in stripped_header}
        date_list = []

        # Read the rest of the lines (data rows)
        for line in file:
            values = line.strip().split(";")
            values = values[:4] + values[5:]
            # Convert the values to floats
            row = {col.strip(): float(value) for col, value in zip(stripped_header, values)}
            if row["hauteur"] != 0.0:
                non_zero_height_encountered = True
                
                # Extract date values
                year = int(values[ian_idx])
                month = int(values[mo_idx])
                day = int(values[jo_idx])
                date_str = f"{year:04d}-{month:02d}-{day:02d}"
                date_list.append(date_str)

                for col in stripped_header:
                    data_dict[col.strip()].append(row[col.strip()])
            if non_zero_height_encountered and (row["hauteur"] == 0.0):
                break

    # start = 21 # 23
    # end = 140
    # density = 10 # density = 20 plants/m2 = 0.002 plants/cm2

    # Thermal time
    thermal_time = [float(i) for i in data_dict["somupvtsem"]]
    # thermal_time = list(np.cumsum([float(i) for i in data_dict["tempeff"]]))
    # thermal_time = list(np.cumsum([float(i) for i in data_dict["tmoy(n)"][:end]]))

    # Green LAI
    plant_leaf_area = [10000*float(i)/density for i in data_dict["laimax"]] # from m2/m2 to cm2/plant
    leaf_area_incr = [plant_leaf_area[0]] + [plant_leaf_area[i+1]-plant_leaf_area[i] for i in range(len(plant_leaf_area[1:]))]

    # Senescent LAI
    sen_leaf_area = [10000*float(i)/density for i in data_dict["laisen(n)"]] # from m2/m2 to cm2/plant
    sen_leaf_area_incr = [sen_leaf_area[0]] + [sen_leaf_area[i+1]-sen_leaf_area[i] for i in range(len(sen_leaf_area[1:]))]

    # Phenology
    emergence = data_dict["ilevs"][-1] - data_dict["jul"][0] # from pseudo julian day (from the beginning of the year) to day from begining of the simulation
    end_juv = data_dict["iamfs"][-1] - data_dict["jul"][0]
    max_lai = data_dict["ilaxs"][-1] - data_dict["jul"][0]

    # Height
    height = [float(i)*100 for i in data_dict["hauteur"]] # from m to cm
    height_incr = [height[0]] + [height[i+1]-height[i] for i in range(len(height[1:]))]

    # Incident PAR
    par_inc = [0.95*0.48*float(j) for j in data_dict["trg(n)"]]
    # Absorbed PAR
    par_abs = [float(i)/(0.95*0.48*float(j)) for i, j in zip(data_dict["raint"], data_dict["trg(n)"])] # to % of light intercepted, in MJ/m^2

    return {
        i+1: {"Date": date_list[i],
            "Thermal time": round(thermal_time[i],4),
            "Phenology": 'germination' if i+1 < emergence else 'juvenile' if emergence <= i+1 < end_juv else 'exponential' if end_juv <= i+1 < max_lai else 'repro',
            "Plant leaf area": round(plant_leaf_area[i],4), 
            "Leaf area increment": round(leaf_area_incr[i],4), 
            "Plant senescent leaf area": round(sen_leaf_area[i],4),
            "Senescent leaf area increment": round(sen_leaf_area_incr[i],4),
            "Plant height": round(height[i],4), 
            "Height increment": round(height_incr[i],4), 
            "Incident PAR": round(par_inc[i],4),
            "Absorbed PAR": round(par_abs[i],4)}
        for i in range(len(thermal_time))
    }


def get_stics_management_params(file_tec_xml):
    """Retrieve STICS management parameters from an XML file."""
    params_tec = ['densitesem', 'interrang']
    return read_xml_file(file_tec_xml, params_tec)

def get_stics_senescence_params(file_plt_xml):
    """Retrieve STICS senescence parameters from an XML file."""
    params_sen = ['durvieF', 'ratiodurvieI']
    return read_xml_file(file_plt_xml, params_sen)

def get_stics_dynamics(stics_output_file, sowing_density):
    """Retrieve STICS growth and senescence dynamics from a STICS output file."""
    return read_sti_file(stics_output_file, sowing_density)

def get_stics_data(file_tec_xml, file_plt_xml, stics_output_file):
    """Retrieve STICS management and senescence parameters, and growth dynamics."""
    tec_stics = get_stics_management_params(file_tec_xml)
    sowing_density = tec_stics['densitesem']
    
    stics_output_data = get_stics_dynamics(stics_output_file, sowing_density)
    
    sen_stics = get_stics_senescence_params(file_plt_xml)
    lifespan = sen_stics['durvieF']
    lifespan_early = sen_stics['ratiodurvieI'] * lifespan
    
    return sowing_density, stics_output_data, lifespan, lifespan_early


def stics_weather_3d(filename, daily_dynamics):
    """Load the weather data from a file and filter it based on the first and last dates of plant growth."""
    df = meteo_day(filename)  # noqa: PD901

    # Get the first and last dates from daily_dynamics
    first_date = list(daily_dynamics.values())[0]["Date"]  # noqa: RUF015
    last_date = list(daily_dynamics.values())[-1]["Date"]

    # Use these dates to filter your DataFrame
    return df[(df.daydate >= pd.to_datetime(first_date)) & (df.daydate <= pd.to_datetime(last_date))]
