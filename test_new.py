#!/usr/bin/env python3

from pathlib import Path
import psutil
from time import sleep
import time
from threading import Thread
from math import sqrt
import pandas as pd
import numpy as np
from osgeo import gdal
import json
from shapely.geometry import Point, Polygon
gdal.UseExceptions()


class ComputeThread(Thread):
    def __init__(self, side, df, L):
        Thread.__init__(self)
        self.side = side
        self.df = df
        self.L = L
        self.df = df
    def get_mem_usg(self):
        return psutil.virtual_memory().percent
    def check_pause(self):
        mem_usg = self.get_mem_usg()
        while mem_usg >= 75:
            mem_usg = self.get_mem_usg()
            print(f"Memory usage is {mem_usg}. Pausing {self.triangle} til <= 75%")
            sleep(30)
    def run(self):
        self.check_pause()
        print(f"Running {self.side}")
        if self.side == "horiz":
            df_out = abs(self.df - self.df.shift(axis=1, periods=-1)).applymap(
                    lambda a: sqrt(a**2 + self.L**2)/2)
        elif self.side == "vert":
            df_out = abs(self.df - self.df.shift(axis=0, periods=-1)).applymap(
                    lambda a: sqrt(a**2 + self.L**2)/2)
        elif self.side == "l_up":
            df_out = abs(self.df - self.df.shift(axis=0, periods=-1).shift(
                axis=1, periods=-1)).applymap(lambda a: sqrt(a**2 + (sqrt(2*self.L**2))**2)/2)
        elif self.side == "r_up":
            df_out = abs(self.df.shift(axis=1, periods=-1) - self.df.shift(axis=0, periods=-1)
                    ).applymap(lambda a: sqrt(a**2 + (sqrt(2*self.L**2))**2)/2)
        else:
            print(f"Pleas choose: horiz, vert, l_up, or r_up. Side supplied was {self.side}")
        df_out.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
        self.df_out = df_out.astype('single')

class SAThread(Thread):
    def __init__(self, triangle, a, b, c):
        Thread.__init__(self)
        self.triangle = triangle
        x, y = a.shape
        self.a = a.iloc[0:x-2, 0:y-2]
        self.b = b.iloc[0:x-2, 0:y-2]
        self.c = c.iloc[0:x-2, 0:y-2]
        self.sa = 0
    def get_mem_usg(self):
        return psutil.virtual_memory().percent
    def check_pause(self):
        mem_usg = self.get_mem_usg()
        while mem_usg >= 75:
            mem_usg = self.get_mem_usg()
            print(f"Memory usage is {mem_usg}. Pausing {self.triangle} til <= 75%")
            sleep(30)
    def run(self):
        self.check_pause()
        s = (self.a + self.b + self.c)/2
        base = s*(s - self.a)*(s - self.b)*(s - self.c)
        sa = base.applymap(lambda x: sqrt(abs(x)))
        sa.replace([np.inf, -np.inf], 0, inplace=True)
        self.sa = round(sa.sum().sum(), 2)
        print(f"Surface area for {self.triangle}: {self.sa}")


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
    a = Point(xs[0], ys[0]).within(poly)
    b = Point(xs[0], ys[-1]).within(poly)
    c = Point(xs[-1], ys[0]).within(poly)
    d = Point(xs[-1], ys[-1]).within(poly)
    if a and b and c and d:
        print("Map fully within the boundary")
        return pd.DataFrame([[True]*ncols]*nrows), True
    elif not a and not b and not c and not d:
        print("Map not in the boundary")
        return pd.DataFrame([[False]*ncols]*nrows), False
    print("Returning a parital map")
    cols = []
    for x in xs:
        row = []
        for y in ys:
            row.append((x,y))
        cols.append(row)
    df = pd.DataFrame(cols)
    return df.applymap(lambda a: Point(a).within(poly)), True


start = time.time()
total_sa = 0
f2 = open("states.json")
states = json.load(f2)
poly = Polygon(states["Idaho"])
f = 'DEM/USGS_1_n46w116_20220309.tif'
f = 'DEM/USGS_1_n42w112_20211101.tif'
directory = "DEM"
files = Path(directory).glob('*')
for f in files:
    print(f"map: {f}")
    ds = gdal.Open(str(f))
    elevation = ds.ReadAsArray()
    coords, ok = get_real_coords(elevation, poly)
    if not ok:
        continue
    #elevation = [[190,170,155],[183,165,145],[175,160,122]]
    #L = 100
    L = 30
    df = pd.DataFrame(elevation)
    # apply coordinate mask
    df = df[coords]
    df = df.astype('float16')

    dataset_sa = 0
    threads = {}
    sides = {"horiz": None, "vert": None, "l_up": None, "r_up": None}
    for side in sides:
        threads[side] = ComputeThread(side, df, L)
        threads[side].start()
    for side in sides:
        threads[side].join()
        sides[side] = threads[side].df_out

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
        threads[name] = SAThread(name, tri["horiz"], tri["vert"], tri["diag"])
        threads[name].start()
        sleep(30)
    for name in names:
        threads[name].join()
        dataset_sa += threads[name].sa


    end = time.time()
    total_sa += dataset_sa
    print(f"file: {f}, dataset_sa: {dataset_sa}, elapsed_time: {round(end-start,0)}")
print(f"total_sa: {total_sa}")

