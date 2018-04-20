#!/usr/bin/env pythonw

# -*- coding: utf-8 -*-

import os
import re
import sys

import numpy as np

# I/O functions ... boring
def pw_char(x):
    r = x % 1600

    if r < 401:
        return '|'
    elif r < 801:
        return '\\'
    elif r < 1201:
        return '—'
    else:
        return '/'

def load_status(i, maximum):
    sys.stdout.write("Loading ({1}%)... ({0})\r".format(pw_char(i), int(np.floor(100 * i / maximum))))

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
            load_status(currentline, linecount)

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

def load_cache(filepath):
    return np.loadtxt(open(filepath + ".cache"), delimiter=",")

def write_cache(filepath, data):
    np.savetxt(filepath + ".cache", data, delimiter=",")

def parse_line(line, dtype=np.float):
    vals = re.split("\s|\,\s", line)

    # Clean up whitespace and return as np array with correct datatype
    return np.array(list(filter(lambda s: s != "", vals))).astype(dtype)

# Physics functions... interesting :D
def get_alpha_kb_ratio(data, temp, column=0, units=1.):
    #    a*<x^2> = kbT
    # -> a / kb  = T / <x^2>

    # Subtract the mean value first (since <x^2> is assumed from x_0 = 0)
    x = data[:, column] - np.mean(data[:, column])
    
    # Return with the corrected variance
    return temp / np.var(units * x)

##### EXECUTION #####
T   = 293.15              # Assumed to be room temperature

convertXdata1 = (1. / 573.26e-3) * 1e-6
convertYdata1 = (1. / 550.22e-3) * 1e-6

convertXdata3 = (1. / 683.89e-3) * 1e-6
convertYdata3 = (1. / 559.76e-3) * 1e-6

fn1 = "../data/data1.dat"

dat1, fromcache1 = read_datfile(fn1, skiprows=2)

akbx1 = get_alpha_kb_ratio(dat1, T, units=convertXdata1)
akby1 = get_alpha_kb_ratio(dat1, T, column=1, units=convertYdata1)

print("(\\alpha / k_B)_x = {0:1.4e}\n(\\alpha / k_B)_y = {1:1.4e}".format(akbx1, akby1)) # Off

fn2 = "../data/data2.dat"

dat2, fromcache2 = read_datfile(fn2, skiprows=2)

akbx2 = get_alpha_kb_ratio(dat2, T, units=convertXdata1)
akby2 = get_alpha_kb_ratio(dat2, T, column=1, units=convertYdata1)

print("(\\alpha / k_B)_x = {0:1.4e}\n(\\alpha / k_B)_y = {1:1.4e}".format(akbx2, akby2)) # WAY off

fn3 = "../data/data3.dat"

dat3, fromcache3 = read_datfile(fn3, skiprows=2)

akbx3 = get_alpha_kb_ratio(dat3, T, units=convertXdata3)
akby3 = get_alpha_kb_ratio(dat3, T, column=1, units=convertYdata3)

print("(\\alpha / k_B)_x = {0:1.4e}\n(\\alpha / k_B)_y = {1:1.4e}".format(akbx3, akby3)) # OK at least these two match in order of magnitude

if not fromcache1:
    write_cache(fn1, dat1)

if not fromcache2:
    write_cache(fn2, dat2)

if not fromcache3:
    write_cache(fn3, dat3)

