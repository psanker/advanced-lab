# -*- coding: utf-8 -*-

import numpy as np

import os
import glob


import sys

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

def load_data(filename):
    # dirpath = os.path.normpath(os.path.join(os.getcwd(), "../data"))
    # oldpath = os.getcwd()

    # Set CWD and search files
    # os.chdir(dirpath)

    files = sorted(glob.glob("../data/{}*".format(filename)))

    for i in range(len(files)):
        files[i] = os.path.normpath(os.path.join(os.getcwd(), files[i]))

    count, times, A, B, AB, C, AC, BC, ABC = parse_files(files, filename)

    # Reset CWD to preserve state
    # os.chdir(oldpath)

    sys.stdout.write("Loading data... Done.\n")
    return count, times, A, B, AB, C, AC, BC, ABC

def parse_files(files, filename):
    FILE_START = True
    PROGRESS  = 0

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
    TYPES = (np.int, "|S10", np.int, np.int, np.int, np.int, np.int, np.int, np.int)

    for f in files:
        # Top of stack
        if FILE_START:
            data = np.loadtxt(f, skiprows=4, delimiter=",", dtype={"names": NAMES, "formats": TYPES})
            
            for i in range(len(data)):
                # Update loading tick
                sys.stdout.write("Loading data... {}\r".format(pw_char(PROGRESS)))
                PROGRESS += 1

                count.append(data[i][0])
                times.append(data[i][1])
                countA.append(data[i][2])
                countB.append(data[i][3])
                countAB.append(data[i][4])
                countC.append(data[i][5])
                countAC.append(data[i][6])
                countBC.append(data[i][7])
                countABC.append(data[i][8])

            FILE_START = False
        else:
            data = np.loadtxt(f, delimiter=",", dtype={"names": NAMES, "formats": TYPES})
            
            for i in range(len(data)):
                # Update loading tick
                sys.stdout.write("Loading data... {}\r".format(pw_char(PROGRESS)))
                PROGRESS += 1

                count.append(data[i][0])
                times.append(data[i][1])
                countA.append(data[i][2])
                countB.append(data[i][3])
                countAB.append(data[i][4])
                countC.append(data[i][5])
                countAC.append(data[i][6])
                countBC.append(data[i][7])
                countABC.append(data[i][8])

    return np.array(count), np.array(times), np.array(countA), np.array(countB), np.array(countAB), np.array(countC), np.array(countAC), np.array(countBC), np.array(countABC)

# EXECUTABLE SECTION

print("First run -----")
FILENAME = "firstrun"

count, times, A, B, AB, C, AC, BC, ABC = load_data(FILENAME)

print(np.sum(A), np.sum(B), np.sum(AB))

print("Second run -----")
FILENAME = "secondrun"

count, times, A, B, AB, C, AC, BC, ABC = load_data(FILENAME)

print(np.sum(A), np.sum(B), np.sum(AB))
