#!/usr/bin/env pythonw

# -*- coding: utf-8 -*-

import numpy as np
import re
import sys

def pw_char(x):
    r = x % 12

    if r == 0 or r == 1 or r == 2:
        return '|'
    elif r == 3 or r == 4 or r == 5:
        return '\\'
    elif r == 6 or r == 7 or r == 8:
        return '—'
    else:
        return '/'

def load_status(i):
    sys.stdout.write("Loading... ({})\r".format(pw_char(i)))

def read_datfile(filepath, skiprows=0, dtype=np.float):
    data        = None # Return value to be determined once number of columns is determined

    # Gross overhead, but has to be done
    linecount   = sum(1 for line in open(filepath))

    with open(filepath) as f:
        print("Reading \"{}\"".format(filepath))

        linecount  -= skiprows
        currentline = 1                 # 1-based indexing

        for line in f:
            load_status(currentline)

            if currentline < skiprows + 1:
                currentline += 1
                continue
            else:
                line_data = parse_line(line, dtype=dtype)
                i         = currentline - skiprows - 1

                if i == 0:
                    # First line in processing
                    data = np.empty((linecount, len(line_data)), dtype=dtype)

                    # Insert line to first entry
                    data[0, :] = line_data
                else:
                    data[i, :] = line_data

                currentline += 1
    return data

def parse_line(line, dtype=np.float):
    vals = re.split("\s|\,\s", line)

    # Clean up whitespace
    vals = np.array(list(filter(lambda s: s != "", vals)))

    return vals.astype(dtype)

##### EXECUTION #####
data = read_datfile("../data/data1.dat", skiprows=2)
print(data)

