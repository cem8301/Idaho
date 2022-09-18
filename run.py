#!/usr/bin/env python3

import requests
from pathlib import Path
import psutil
from time import sleep
import time
from math import sqrt
import pandas as pd
import numpy as np
from osgeo import gdal
import json
from shapely.geometry import Point, Polygon
import multiprocessing
import wget
import os
gdal.UseExceptions()


def get_mem_usg():
    return psutil.virtual_memory().percent

def check_pause(item):
    mem_usg = get_mem_usg()
    while mem_usg >= 75:
        mem_usg = get_mem_usg()
        print(f"Memory usage is {mem_usg}. Pausing {item} til <= 75%")
        sleep(30)

def sides_process(side, df, L_vert, L_horiz):
    check_pause(side)
    print(f"Running {side}")
    # Where L_lats
    if side == "horiz":
        df_out = abs(df - df.shift(axis=1, periods=-1)).applymap(
                lambda a: sqrt(a**2 + L_horiz**2)/2)
    elif side == "vert":
        df_out = abs(df - df.shift(axis=0, periods=-1)).applymap(
                lambda a: sqrt(a**2 + L_vert**2)/2)
    elif side == "l_up":
        df_out = abs(df - df.shift(axis=0, periods=-1).shift(axis=1, periods=-1)).applymap(
                lambda a: sqrt(a**2 + sqrt(L_vert**2 + L_horiz**2)**2)/2)
    elif side == "r_up":
        df_out = abs(df.shift(axis=1, periods=-1) - df.shift(axis=0, periods=-1)).applymap(
                lambda a: sqrt(a**2 + (sqrt(L_vert**2 + L_horiz**2))**2)/2)
    else:
        print(f"Please choose: horiz, vert, l_up, or r_up. Side supplied was {side}")
    df_out.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
    return_dict[f"sides_{side}"] = df_out.astype('single')

def sa_process(triangle, a, b, c):
    check_pause(triangle)
    s = (a + b + c)/2
    base = s*(s - a)*(s - b)*(s - c)
    sa = base.applymap(lambda x: sqrt(abs(x)))
    sa.replace([np.inf, -np.inf], 0, inplace=True)
    sa = round(sa.sum().sum(), 2)
    print(f"Surface area for {triangle}: {sa}")
    return_dict[f"sa_{triangle}"] = sa


def coords_process(df, x, cores_range):
    check_pause(x)
    cols = df.iloc[:, cores_range[x]:cores_range[x+1]]
    return_dict[f"cores_{x}"] = cols.applymap(lambda a: Point(a).within(poly))


def get_real_coords(elevation, poly):
    # Get location data:
    nrows, ncols = elevation.shape
    # I'm making the assumption that the image isn't rotated/skewed/etc.
    # This is not the correct method in general, but let's ignore that for now
    # If dxdy or dydx aren't 0, then this will be incorrect
    x0, dx, dxdy, y0, dydx, dy = ds.GetGeoTransform()
    # First check if all corners are in the map. if so, mark whole box as true
    xs = np.arange(0, ncols) * dx + x0
    ys = np.arange(0, nrows) * dy + y0
    x_mid = int(len(xs)/2)
    y_mid = int(len(ys)/2)
    lat_mid = ys[y_mid]
    a = Point(xs[0], ys[0]).within(poly)
    b = Point(xs[0], ys[-1]).within(poly)
    c = Point(xs[-1], ys[0]).within(poly)
    d = Point(xs[-1], ys[-1]).within(poly)
    e = Point(xs[x_mid], ys[-1]).within(poly)
    f = Point(xs[x_mid], ys[0]).within(poly)
    g = Point(xs[-1], ys[y_mid]).within(poly)
    h = Point(xs[0], ys[y_mid]).within(poly)
    if a and b and c and d and e and f and g and h:
        print("Map fully within the boundary")
        return pd.DataFrame([[True]*ncols]*nrows), lat_mid, True
    elif not a and not b and not c and not d and not e and not f and not g and not h:
        print("Map not in the boundary")
        return pd.DataFrame([[False]*ncols]*nrows), lat_mid, False
    print("Returning a parital map")
    cols = []
    for x in xs:
        row = []
        for y in ys:
            row.append((x,y))
        cols.append(row)
    df = pd.DataFrame(cols)
    #df_map = df.applymap(lambda a: Point(a).within(poly))
    ps = {}
    # find greatest multiple
    for multiple in reversed(range(num_cores)):
        if len(xs) % multiple == 0:
            break
    cores_range = np.arange(0, multiple+1) * int(len(xs)/(multiple))
    for x in range(0, multiple):
        ps[x] = multiprocessing.Process(target=coords_process, args=(df,x,cores_range,))
        ps[x].start()
        sleep(1)
    for x in range(0, multiple):
        ps[x].join()
        df.iloc[:, cores_range[x]:cores_range[x+1]] = return_dict[f"cores_{x}"]
    return df, lat_mid, True


start = time.time()
num_cores = psutil.cpu_count()-1
manager = multiprocessing.Manager()
return_dict = manager.dict()
m2_to_mile2 = 3.86102e-7
f2 = open("states.json")
states = json.load(f2)
total_sa = 0
total_flat_sa = 0
for state in states:
    state_sa = 0
    #state = "Idaho"
    print(f"Finding surface area within {state} boundary...............................")
    poly = Polygon(states[state])
    bounds = poly.bounds
    # ref: sco.wisc.edu/2022/01/21/how-big-is-a-degree/
    lat_miles = 69.4
    long_miles = np.cos(np.mean([bounds[i] for i in [1,3]]) * np.pi/ 180) * lat_miles #69.172
    flat_sa = poly.area * lat_miles * long_miles 
    r = requests.get(
        f"https://tnmaccess.nationalmap.gov/api/v1/products?datasets=National%20Elevation"
        "%20Dataset%20(NED)%201%20arc-second&dateType=lastUpdated&bbox="
        f"{bounds[0]},{bounds[1]},{bounds[2]},{bounds[3]}&max=1000"
    )
    boundaries = r.json()
    directory = "DEM"
    # {region: date}
    files_state = {
        each_path.parts[1].split("_")[2]: int(each_path.parts[1].split("_")[3].split(".")[0])
        for each_path in Path(directory).glob('*')
    }
    regions_done = []
    boundary_items = boundaries.get("items", [])
    print(f"boundary items: {len(boundary_items)}")
    for fnum, boundary in enumerate(boundary_items):
        fname = boundary["downloadURL"].split("/")[-1]
        cur_region = fname.split("_")[2]
        cur_date = int(fname.split("_")[3].split(".")[0])
        print(f"Checking if {cur_region} is in {regions_done}")
        if cur_region in regions_done:
            print(f"{cur_region} already calculated. Continue")
            continue
        if Path(f"{directory}/{fname}").is_file():
            f = f"{directory}/{fname}"
            print(f"File already local: {f}")
        else:
            if cur_region in files_state.keys():
                orig_date = files_state[cur_region]
                if orig_date < cur_date:
                    os.remove(f"{directory}/USGS_1_{cur_region}_{orig_date}.tif")
                    file_state[cur_region] = cur_date
                    f = wget.download(boundary["downloadURL"], out=directory)
                    print(f"Found a newer file {orig_date} < {cur_date}. wget: {f}")
                else:
                    f = f"{directory}/USGS_1_{cur_region}_{orig_date}.tif"
                    print(f"Current file is the latest {orig_date} > {cur_date}: {f}")
            else:
                f = wget.download(boundary["downloadURL"], out=directory)
                print(f"wget new file: {f}")
        regions_done.append(cur_region)
        ds = gdal.Open(str(f))
        elevation = ds.ReadAsArray()
        coords, lat_mid, ok = get_real_coords(elevation, poly)
        if not ok:
            continue
        #elevation = [[190,170,155],[183,165,145],[175,160,122]]
        #L = 100
        # Same calculation as above for lat and long. Where vertical lengths remain the same
        # (between lattitudes) and hoirzotnal lines (between longitudes) change
        L_vert = 30.87
        L_horiz = np.cos(lat_mid * np.pi/ 180) * L_vert
        boundary_flat_sa = L_vert * L_horiz * len(elevation)**2
        df = pd.DataFrame(elevation)
        # apply coordinate mask
        df = df[coords]
        df = df.astype('float16')
        dataset_sa = 0
        threads = {}
        sides = {"horiz": None, "vert": None, "l_up": None, "r_up": None}
        for side in sides:
            threads[side] = multiprocessing.Process(
                target=sides_process, args=(side, df, L_vert, L_horiz)
            )
            threads[side].start()
        for side in sides:
            threads[side].join()
            sides[side] = return_dict[f"sides_{side}"]

        # triangle i: AB, BE, EA (l_up)
        # triangle ii: BC, BE, EC (r_up)
        # triangle iii: DE, AD, EA (l_up)
        # triangle iv: EF, CF, EC (r_up)
        # triangle v: DE, DG, EG (r_up)
        # triangle vi: EF, FI, EI (l_up)
        # triangle vii: GH, EH, EG (r_up)
        # triangle viii: HI, EH, EI (l_up)
        tris = [
            {"name": "i",
                "horiz": sides["horiz"],
                "vert": sides["vert"].shift(axis=1, periods=-1),
                "diag": sides["l_up"]},
            {"name": "ii",
                "horiz": sides["horiz"].shift(axis=1, periods=-1),
                "vert": sides["vert"].shift(axis=1, periods=-1),
                "diag": sides["r_up"].shift(axis=1, periods=-1)},
            {"name": "iii",
                "horiz": sides["horiz"].shift(axis=0, periods=-1),
                "vert": sides["vert"],
                "diag": sides["l_up"]},
            {"name": "iv",
                "horiz": sides["horiz"].shift(axis=1, periods=-1).shift(axis=0, periods=-1),
                "vert": sides["vert"].shift(axis=1, periods=-2),
                "diag": sides["r_up"].shift(axis=1, periods=-1)},
            {"name": "v",
                "horiz": sides["horiz"].shift(axis=0, periods=-1), 
                "vert": sides["vert"].shift(axis=0, periods=-1),
                "diag": sides["r_up"].shift(axis=0, periods=-1)},
            {"name": "vi",
                "horiz": sides["horiz"].shift(axis=1, periods=-1).shift(axis=0, periods=-1),
                "vert": sides["vert"].shift(axis=1, periods=-2).shift(axis=0, periods=-1),
                "diag": sides["l_up"].shift(axis=1, periods=-1).shift(axis=0, periods=-1)},
            {"name": "vii",
                "horiz": sides["horiz"].shift(axis=1, periods=-1),
                "vert": sides["vert"].shift(axis=1, periods=-1).shift(axis=0, periods=-1),
                "diag": sides["r_up"].shift(axis=0, periods=-1)},
            {"name": "viii",
                "horiz": sides["horiz"].shift(axis=1, periods=-1).shift(axis=0, periods=-2),
                "vert": sides["vert"].shift(axis=1, periods=-1).shift(axis=0, periods=-1),
                "diag": sides["l_up"].shift(axis=1, periods=-1).shift(axis=0, periods=-1)}
        ]

        threads = {}
        names = []
        for idx, tri in enumerate(tris):
            name = tri["name"]
            names.append(name)
            threads[name] = multiprocessing.Process(
                target=sa_process,
                args=(name, tri["horiz"], tri["vert"], tri["diag"])
            )
            threads[name].start()
            sleep(1)
        for name in names:
            threads[name].join()
            dataset_sa += return_dict[f"sa_{name}"] 


        end = time.time()
        state_sa += dataset_sa
        print(
            f"file ({fnum}/{len(boundaries['items'])} {state}): {f}, "
            f"dataset_sa: {round(dataset_sa)}/{round(boundary_flat_sa)} m^2 "
            f"({round(dataset_sa * m2_to_mile2)}/{round(boundary_flat_sa * m2_to_mile2)} miles^2), "
            f"elapsed_time: {round(end-start,0)}"
        )
    total_sa += state_sa
    total_flat_sa += flat_sa
    print(
        f"state_sa {state}: {round(state_sa)} m^2 ({round(state_sa * m2_to_mile2)} miles^2), "
        f"flat_sa_calc: {round(flat_sa / m2_to_mile2)} m^2 ({round(flat_sa)}) miles^2" 
    )
print(
    f"total_sa lower 48: {round(total_sa)} m^2 ({round(total_sa * m2_to_mile2)} miles^2), "
    f"flat_sa_calc: {round(total_flat_sa / m2_to_mile2)} m^2 ({round(total_flat_sa)}) miles^2" 
)
