#!/usr/bin/env pythonw

# -*- coding: utf-8 -*-

#### 1. IMPORTS
import numpy as np

import os
import glob
from datetime import datetime
from astropy.table import Table, Column

import sys

#### 2. FUNCTIONS
# Little loading indicator
def pw_char(i):
    if i % 4 == 0:
        return "(-)"
    elif i % 4 == 1:
        return "(/)"
    elif i % 4 == 2:
        return "(|)"
    else:
        return "(\\)"

# In all honesty, using Pandas would be better for handling this. However..
def load_data(filename):
    files = sorted(glob.glob("../data/{}*".format(filename)))

    for i in range(len(files)):
        files[i] = os.path.normpath(os.path.join(os.getcwd(), files[i]))

    times, data = parse_files(files, filename)

    sys.stdout.write("Loading data... Done.\n")
    return times, data

def parse_files(files, filename):
    FILE_START = True

    count    = []
    times    = []
    countA   = []
    countB   = []
    countAB  = []
    countC   = []
    countAC  = []
    countBC  = []
    countABC = []

    NAMES = ("Count", "Time", "A", "B", "AB", "C (Noise)", "AC", "BC", "ABC")
    TYPES = (np.int, "|S10", np.float, np.float, np.float, np.float, np.float, np.float, np.float)

    for f in files:
        # Top of stack
        data = None

        if FILE_START:
            data = np.loadtxt(f, skiprows=4, delimiter=",", dtype={"names": NAMES, "formats": TYPES})
            FILE_START = False
        else:
            data = np.loadtxt(f, delimiter=",", dtype={"names": NAMES, "formats": TYPES})

        if data is None:
            raise Exception("Data not loaded properly")

        for i in range(len(data)):
            # Update loading tick
            sys.stdout.write("Loading data... {}\r".format(pw_char(i)))

            count.append(data[i][0])
            times.append(data[i][1])
            countA.append(data[i][2])
            countB.append(data[i][3])
            countAB.append(data[i][4])
            countC.append(data[i][5])
            countAC.append(data[i][6])
            countBC.append(data[i][7])
            countABC.append(data[i][8])

    return np.array(times), np.transpose(np.array([countA, countB, countAB, countC, countAC, countBC, countABC]))

def runtime(times):
    start = times[0].decode("utf-8")
    end   = times[-1].decode("utf-8")

    start = datetime.strptime(start, "%I:%M:%S %p")
    end   = datetime.strptime(end, "%I:%M:%S %p")

    return (end - start).total_seconds()

def g2o(times, data, counttime=(1. / 6.554e6)):
    # N_AB / (N_A * N_B)
    count_frac = (np.sum(data[:, 2]) - np.sum(data[:, 6])) / (np.sum(data[:, 0]) * np.sum(data[:, 1]))

    # Total time / read time
    time_frac = runtime(times) / counttime # From settings

    return count_frac * time_frac

def get_counts(times, data):
    return np.array([np.sum(data[:, 0]), np.sum(data[:, 1]), np.sum(data[:, 2]), np.sum(data[:, 3]), np.sum(data[:, 4]), np.sum(data[:, 5]), np.sum(data[:, 6])])

#### 3. EXECUTION

print("First run -----")
FILENAME = "firstrun"

t1, d1 = load_data(FILENAME)

g2o1 = g2o(t1, d1)
counts1 = get_counts(t1, d1)

print("Second run -----")
FILENAME = "secondrun"

t2, d2 = load_data(FILENAME)

g2o2 = g2o(t2, d2)
counts2 = get_counts(t2, d2)

arr = np.array([g2o1, g2o2])
counts = np.array([counts1, counts2]).T
times = np.array([runtime(t1), runtime(t2)]).T

print("\ng2(0): {0:1.4e} \pm {1:1.4e}\n".format(np.mean(arr), np.std(arr)))
print(Table({
    'Counts A': counts[0],
    'Counts B': counts[1],
    'Counts AB': counts[2],
    'Counts C': counts[3],
    'Counts AC': counts[4],
    'Counts BC': counts[5],
    'Counts ABC': counts[6],
    'Runtime': times[:] / 3600.0 #(seconds to hours)
}))
