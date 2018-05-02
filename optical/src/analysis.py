#!/usr/bin/env pythonw
# -*- coding: utf-8 -*-

import os
import re

from concurrent.futures import ProcessPoolExecutor, wait

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy import units as u
from scipy.optimize import curve_fit

# I/O functions ... boring
def read_datfile(filepath, skiprows=0, dtype=np.float64, skipcache=False):
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
def propagate(f, arrx, arrsx, dx=1e-5):
    assert len(arrx) == len(arrsx)

    dfdx2 = np.zeros((len(arrx), len(arrx)))

    for i in range(len(arrx)):
        temp = arrx
        temp[i] += dx*1j

        dfdx2[i, i] = np.imag(f(temp) / dx)**2

    return (arrsx.T) @ dfdx2 @ arrsx

def get_alpha_kb_ratio(data, temp, column=0, units=1., sunits=1., sumcolumn=2):
    #    a*<x^2> = kbT
    # -> a / kb  = T / <x^2>

    # Subtract the mean value first (since <x^2> is assumed from x_0 = 0)
    x = data[:, column] - np.mean(data[:, column])

    # According to manual, X and Y may be normalized by SUM. So, de-normalize; also, add units
    x *= data[:, sumcolumn] * u.V

    # Corrected variance
    mu = temp / (units * units * np.mean(x * x))

    # Return with the corrected variance
    return mu

def fit_power_spectrum(data, temp, column=1, sumcolumn=3): # column=1: x power spec; column=2: y power spec
    # Data in the power spectrum data needs to be squared
    # Reject first row to remove average
    pspec = (data[1:, column] * data[1:, sumcolumn])

    #              kbT       <-- p0
    # P^2 = ----------------
    #         p1*f^2 + p2
    def func(f, *p):
        return p[0] / (p[1]*f**2 + p[2])

    p    = np.array([0.3, 0.3, 0.3]) # Bad initial guess... Maybe should provide better guess

    return curve_fit(func, data[1:, 0], pspec, p0=p) # returns (params, cov matrix)

def plot_power_spectrum(data, pX, pY):
    fig, ax = plt.subplots(2)
    ax[0].plot(data[1:, 0], data[1:, 1]*data[1:, 3])

    def func(f, p):
        return p[0] / (p[1]*f**2 + p[2])

    freq = np.linspace(0, np.amax(data[:, 0]), 80000)
    ax[0].plot(freq, func(freq, pX))

    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    
    ax[1].plot(data[1:, 0], data[1:, 2])
    ax[1].plot(freq, func(freq, pY))

    ax[1].set_xscale("log")
    ax[1].set_yscale("log")

def process_frequency_data(filename, opts={}):
    outstring = ""

    data, fromcache = read_datfile(filename, opts["skiprows"])

    paramsX, covX = fit_power_spectrum(data, opts["T"], column=1)
    paramsY, covY = fit_power_spectrum(data, opts["T"], column=2)

    outstring += "\n--------------------\n"
    outstring += "Fourier dataset \"{}\"\n".format(opts["dataname"])

    outstring += "p0_x = {0:1.3e} \\pm {1:1.3e}\n".format(paramsX[0], covX[0, 0]**0.5)
    outstring += "p0_y = {0:1.3e} \\pm {1:1.3e}\n".format(paramsY[0], covY[0, 0]**0.5)
    outstring += "p1_x = {0:1.3e} \\pm {1:1.3e}\n".format(paramsX[1], covX[1, 1]**0.5)
    outstring += "p1_y = {0:1.3e} \\pm {1:1.3e}\n".format(paramsY[1], covY[1, 1]**0.5)
    outstring += "p2_x = {0:1.3e} \\pm {1:1.3e}\n".format(paramsX[2], covX[2, 2]**0.5)
    outstring += "p2_y = {0:1.3e} \\pm {1:1.3e}".format(paramsY[2], covY[2, 2]**0.5)

    # This is broken lul
    # if not fromcache:
    #     write_cache(filename, data)

    return outstring, data, paramsX, paramsY

def process_position_data(filename, opts={}):
    outstring = ""

    data, fromcache = read_datfile(filename, opts["skiprows"])

    akbx = get_alpha_kb_ratio(data, opts["T"], units=opts["convertX"])
    akby = get_alpha_kb_ratio(data, opts["T"], column=1, units=opts["convertY"])

    outstring += "\n--------------------\n"
    outstring += "Position dataset \"{}\"\n".format(opts["dataname"])
    outstring += "(\\alpha / k_B)_x = {0:1.4e}\n(\\alpha / k_B)_y = {1:1.4e}".format(akbx, akby)

    if not fromcache:
        write_cache(filename, data)

    return outstring

##### EXECUTION #####
T = 293.15 * u.K           # Assumed to be room temperature

convertXdata1 = ((1. / 573.26e-3) * (u.micron / u.V)).to(u.m / u.V)
convertYdata1 = ((1. / 550.22e-3) * (u.micron / u.V)).to(u.m / u.V)

convertXdata3 = ((1. / 683.89e-3) * (u.micron / u.V)).to(u.m / u.V)
convertYdata3 = ((1. / 559.76e-3) * (u.micron / u.V)).to(u.m / u.V)

convertXdata4 = ((1. / 745.20e-3) * (u.micron / u.V)).to(u.m / u.V)
convertYdata4 = ((1. / 703.16e-3) * (u.micron / u.V)).to(u.m / u.V)

convertXdata5 = ((1. / 1.09) * (u.micron / u.V)).to(u.m / u.V)
convertYdata5 = ((1. / 1.00) * (u.micron / u.V)).to(u.m / u.V)

def main():
    dataset = []
    pX      = []
    pY      = []

    # By default, the pool executor only spins off 5 threads. This should be enough for us.
    with ProcessPoolExecutor(max_workers=5) as pool:
        futures = []

        futures.append(pool.submit(process_position_data, "../data/data1.dat", opts={
            "dataname": "Data 1",
            "skiprows": 2,
            "convertX": convertXdata1,
            "convertY": convertYdata1,
            "T": T
        }))

        futures.append(pool.submit(process_position_data, "../data/data2.dat", opts={
            "dataname": "Data 2",
            "skiprows": 2,
            "convertX": convertXdata1,
            "convertY": convertYdata1,
            "T": T
        }))

        futures.append(pool.submit(process_position_data, "../data/data3.dat", opts={
            "dataname": "Data 3",
            "skiprows": 2,
            "convertX": convertXdata3,
            "convertY": convertYdata3,
            "T": T
        }))

        futures.append(pool.submit(process_position_data, "../data/data4.dat", opts={
            "dataname": "Data 4",
            "skiprows": 2,
            "convertX": convertXdata4,
            "convertY": convertYdata4,
            "T": T
        }))

        futures.append(pool.submit(process_position_data, "../data/data5.dat", opts={
            "dataname": "Data 5",
            "skiprows": 2,
            "convertX": convertXdata5,
            "convertY": convertYdata5,
            "T": T
        }))

        futures.append(pool.submit(process_frequency_data, "../data/freq1.FDdat", opts={
            "dataname": "Data 1",
            "skiprows": 4,
            "T": T
        }))

        futures.append(pool.submit(process_frequency_data, "../data/freq3.FDdat", opts={
            "dataname": "Data 3",
            "skiprows": 4,
            "T": T
        }))

        futures.append(pool.submit(process_frequency_data, "../data/freq4.FDdat", opts={
            "dataname": "Data 4",
            "skiprows": 4,
            "T": T
        }))

        futures.append(pool.submit(process_frequency_data, "../data/freq5.FDdat", opts={
            "dataname": "Data 5",
            "skiprows": 4,
            "T": T
        }))

        for x in wait(futures)[0]:
            res = x.result()

            # Check if future is from freq space
            if type(res) is tuple:
                print(res[0])
                dataset.append(res[1])
                pX.append(res[2])
                pY.append(res[3])
            else:
                print(res)
    
    plot_power_spectrum(dataset[0], pX[0], pY[0])
    plot_power_spectrum(dataset[1], pX[1], pY[1])
    plot_power_spectrum(dataset[2], pX[2], pY[2])
    plt.show()

if __name__ == "__main__":
    main()
