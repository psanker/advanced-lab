#!/usr/bin/env pythonw

# -*- coding: utf-8 -*-

import os
import re

from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

from astropy import units as u

# I/O functions ... boring
def read_datfile(filepath, skiprows=0, dtype=np.float, skipcache=False):
    data        = None # Return value to be determined once number of columns is determined

    # Determine if a file has been loaded before and therefore has been formatted properly
    if has_cache(filepath) and not skipcache:
        print("Loading \"{}\" from cache...".format(filepath))
        data = load_cache(filepath)
        return data, True # Immediately return

    # Gross overhead, but has to be done in order to get the memory allocation right
    linecount   = sum(1 for line in open(filepath))

    with open(filepath) as f:
        print("Reading \"{}\"".format(filepath))

        linecount  -= skiprows  # The actual maximum number of entries
        currentline = 1         # 1-based indexing

        for line in f:
            if currentline < skiprows + 1:
                currentline += 1
                continue
            else:
                line_data = parse_line(line, dtype=dtype)
                i         = currentline - skiprows - 1

                if i == 0:
                    # First line in processing
                    data = np.empty((linecount, len(line_data)), dtype=dtype)

                data[i, :] = line_data
                currentline += 1
    return data, False

def has_cache(filepath):
    return os.path.exists(filepath + ".cache") and os.path.isfile(filepath + ".cache")

def load_cache(filepath, columns=3):
    return pd.read_csv(filepath + ".cache", delimiter=",").values # SO much faster

def write_cache(filepath, data):
    np.savetxt(filepath + ".cache", data, delimiter=",")

def parse_line(line, dtype=np.float):
    vals = re.split("\s|\,\s", line)

    # Clean up whitespace and return as np array with correct datatype
    return np.array(list(filter(lambda s: s != "", vals))).astype(dtype)

# Physics functions... interesting :D
def get_alpha_kb_ratio(data, temp, column=0, units=1., sumcolumn=2):
    #    a*<x^2> = kbT
    # -> a / kb  = T / <x^2>

    # Subtract the mean value first (since <x^2> is assumed from x_0 = 0)
    x = data[:, column] - np.mean(data[:, column])

    # According to manual, X and Y may be normalized by SUM. So, de-normalize; also, add units
    x *= data[:, sumcolumn] * u.V

    # Return with the corrected variance
    return temp / np.var(units * x)

# Each individual future is an instance of this function
def process_data(filename, opts={}):
    outstring = ""

    data, fromcache = read_datfile(filename, opts["skiprows"])

    akbx = get_alpha_kb_ratio(data, opts["T"], units=opts["convertX"])
    akby = get_alpha_kb_ratio(data, opts["T"], column=1, units=opts["convertY"])

    outstring += "Dataset \"{}\"\n".format(opts["dataname"])
    outstring += "(\\alpha / k_B)_x = {0:1.4e}\n(\\alpha / k_B)_y = {1:1.4e}\n".format(akbx, akby)

    if not fromcache:
        write_cache(filename, data)

    return outstring

##### EXECUTION #####
T   = 293.15 * u.K           # Assumed to be room temperature

convertXdata1 = ((1. / 573.26e-3) * (u.micron / u.V)).to(u.m / u.V)
convertYdata1 = ((1. / 550.22e-3) * (u.micron / u.V)).to(u.m / u.V)

convertXdata3 = ((1. / 683.89e-3) * (u.micron / u.V)).to(u.m / u.V)
convertYdata3 = ((1. / 559.76e-3) * (u.micron / u.V)).to(u.m / u.V)

convertXdata4 = ((1. / 745.20e-3) * (u.micron / u.V)).to(u.m / u.V)
convertYdata4 = ((1. / 703.16e-3) * (u.micron / u.V)).to(u.m / u.V)

convertXdata5 = ((1. / 1.09) * (u.micron / u.V)).to(u.m / u.V)
convertYdata5 = ((1. / 1.00) * (u.micron / u.V)).to(u.m / u.V)

pool = ThreadPoolExecutor(5)

futures = []

futures.append(pool.submit(process_data, "../data/data1.dat", opts={
    "dataname": "Data 1",
    "skiprows": 2,
    "convertX": convertXdata1,
    "convertY": convertYdata1,
    "T": T
}))

futures.append(pool.submit(process_data, "../data/data2.dat", opts={
    "dataname": "Data 2",
    "skiprows": 2,
    "convertX": convertXdata1,
    "convertY": convertYdata1,
    "T": T
}))

futures.append(pool.submit(process_data, "../data/data3.dat", opts={
    "dataname": "Data 3",
    "skiprows": 2,
    "convertX": convertXdata3,
    "convertY": convertYdata3,
    "T": T
}))

futures.append(pool.submit(process_data, "../data/data4.dat", opts={
    "dataname": "Data 4",
    "skiprows": 2,
    "convertX": convertXdata4,
    "convertY": convertYdata4,
    "T": T
}))

futures.append(pool.submit(process_data, "../data/data5.dat", opts={
    "dataname": "Data 5",
    "skiprows": 2,
    "convertX": convertXdata5,
    "convertY": convertYdata5,
    "T": T
}))

for x in as_completed(futures):
    print(x.result())
