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
from scipy.stats import norm

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
def propagate(f, arrx, arrsx, dx=1e-8):
    assert len(arrx) == len(arrsx)

    dfdx2 = np.zeros((len(arrx), len(arrx)))

    for i in range(len(arrx)):
        # Complex step differentiation
        temp = arrx + 0j # Promote array to complex values
        temp[i] += dx*1j # Step the differentiating variable

        dfdx2[i, i] = np.imag(f(temp) / dx)**2 # Evaluate

    return np.matmul(np.matmul((arrsx.T), dfdx2), arrsx)

def get_alpha_kb_ratio(data, temp, column=0, units=1., sunits=1., sumcolumn=2):
    #    a*<x^2> = kbT
    # -> a / kb  = T / <x^2>

    # Subtract the mean value first (since <x^2> is assumed from x_0 = 0)
    x = data[:, column] - np.mean(data[:, column])

    # According to manual, X and Y MUST be normalized by SUM.
    x /= data[:, sumcolumn]

    # Corrected variance
    mu  = temp / (np.mean(x * x)*(u.V)**2 / (units**2))
    smu = np.sqrt(propagate(lambda a: temp.value / (a[0]**(-2) * np.mean(x * x)), np.array([units.value]), np.array([sunits.value])))

    # Return with the corrected variance
    return mu, smu * mu.unit

def fit_power_spectrum(data, temp, column=1, sumcolumn=3): # column=1: x power spec; column=2: y power spec
    # Data in the power spectrum data needs to be squared
    # Reject first row to remove average
    pspec = (data[1:, column])

    #              kbT       <-- p0: Turns out to be not a good parameter to use
    # P^2 = ----------------
    #       p1*(f^2 + p2^2)  <-- p2: The corner frequency is desired
    def func(f, *p):
        return p[0] / (p[1]*(f**2 + p[2]**2))

    p    = np.array([1e-6, 0.3, 1]) # Bad initial guess... Maybe should provide better guess

    return curve_fit(func, data[1:, 0], pspec, p0=p, method="lm") # returns (params, cov matrix)

def interpret_alpha(p2, dp2, a=3.*u.micron):
    # Extract the spring constant from the cornering frequency
    TWO_PI    = 2.*np.pi
    eta       = 8.90e-4 * (u.Pa * u.s)
    beta = 3. * TWO_PI * eta * a

    P2  = p2 * (1 / u.s)
    DP2 = dp2 * (1 / u.s)

    alpha  = TWO_PI * beta * np.abs(P2)
    salpha = np.sqrt(propagate(lambda a: (TWO_PI * beta).value * a[0], np.array([P2.value]), np.array([DP2.value])))

    return alpha.to("N / m"), (salpha * alpha.unit).to("N / m")

def extract_kb_alpha(postup, freqtup):
    akbx, sakbx = postup[1]
    akby, sakby = postup[2]
    alpha_x, salpha_x = freqtup[4]
    alpha_y, salpha_y = freqtup[5]

    kb_x = ((1 / akbx) * alpha_x)
    kb_y = ((1 / akby) * alpha_y)

    skb_x = np.sqrt(propagate(lambda a: (1 / a[0])*a[1], np.array([akbx.value, alpha_x.value]), np.array([sakbx.value, salpha_x.value])))
    skb_y = np.sqrt(propagate(lambda a: (1 / a[0])*a[1], np.array([akby.value, alpha_y.value]), np.array([sakby.value, salpha_y.value])))

    skb_x = (skb_x * kb_x.unit).to("J / K")
    kb_x  = (kb_x).to("J / K")
    skb_y = (skb_y * kb_y.unit).to("J / K")
    kb_y  = (kb_y).to("J / K")

    print("\"{0}\" kB_x: {1:1.3e} \\pm {2:1.3e}".format(postup[0], kb_x, skb_x))
    print("\"{0}\" kB_y: {1:1.3e} \\pm {2:1.3e}".format(postup[0], kb_y, skb_y))

    return kb_x, kb_y, skb_x, skb_y, alpha_x, alpha_y, salpha_x, salpha_y

def plot_kb(exp,unc,theory=1.381e-23):
    x = np.linspace(exp - 10*unc, exp + 10*unc, 100)
    plt.figure()
    plt.plot(x,norm.pdf(x, exp, unc), label=('$k_B$={0:1.3e} $\pm$ {1:1.3e} $J/K$'.format(exp, unc)))
    plt.axvline(theory, ls='--', color='k', label=('$k_B$={0:1.3e} $J/K$'.format(theory)))
    plt.legend(loc='lower left')

def plot_current_dependence(alphas, salphas, currents):
    # assert len(alphas) == 2*len(currents), 'Array dimensions do not match'
    plt.figure()
    plt.errorbar(currents,alphas[0],yerr=salphas[0],fmt='o',markersize=3,label='x')
    plt.errorbar(currents,alphas[1],yerr=salphas[1],fmt='o',markersize=3,label='y')
    plt.ylim(0,1e-5)
    plt.legend()

def plot_power_spectrum(datatuple):
    dname = datatuple[0]
    data  = datatuple[1]
    pX    = datatuple[2]
    pY    = datatuple[3]

    fig, ax = plt.subplots(2)
    ax[0].plot(data[1:, 0], data[1:, 1])

    def func(f, p):
        return p[0] / (p[1]*(f**2 + p[2]**2))

    freq = np.linspace(0, np.amax(data[:, 0]), 80000)
    ax[0].plot(freq, func(freq, pX))

    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("{0} - Frequency ($Hz$)".format(dname))
    ax[0].set_ylabel("Power Spectrum in X ($V^2 / Hz$)")

    ax[1].plot(data[1:, 0], data[1:, 2])
    ax[1].plot(freq, func(freq, pY))

    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("{0} - Frequency ($Hz$)".format(dname))
    ax[1].set_ylabel("Power Spectrum in Y ($V^2 / Hz$)")

def process_frequency_data(filename, opts={}):
    outstring = ""

    data, fromcache = read_datfile(filename, opts["skiprows"])

    paramsX, covX = fit_power_spectrum(data, opts["T"], column=1)
    paramsY, covY = fit_power_spectrum(data, opts["T"], column=2)

    outstring += "\n--------------------\n"
    outstring += "Fourier dataset \"{}\"\n".format(opts["dataname"])

    outstring += "p0_x = {0:1.3e} \\pm {1:1.3e}\n".format(paramsX[0], np.abs(covX[0, 0])**0.5)
    outstring += "p0_y = {0:1.3e} \\pm {1:1.3e}\n".format(paramsY[0], np.abs(covY[0, 0])**0.5)
    outstring += "p1_x = {0:1.3e} \\pm {1:1.3e}\n".format(paramsX[1], np.abs(covX[1, 1])**0.5)
    outstring += "p1_y = {0:1.3e} \\pm {1:1.3e}\n".format(paramsY[1], np.abs(covY[1, 1])**0.5)
    outstring += "p2_x = {0:1.3e} \\pm {1:1.3e}\n".format(paramsX[2], np.abs(covX[2, 2])**0.5)
    outstring += "p2_y = {0:1.3e} \\pm {1:1.3e}\n".format(paramsY[2], np.abs(covY[2, 2])**0.5)

    alpha_x, salpha_x = interpret_alpha(paramsX[2], np.abs(covX[2, 2])**0.5)
    alpha_y, salpha_y = interpret_alpha(paramsY[2], np.abs(covY[2, 2])**0.5)

    outstring += "\\alpha_x = {0:1.3e} \\pm {1:1.3e}\n".format(alpha_x, salpha_x)
    outstring += "\\alpha_y = {0:1.3e} \\pm {1:1.3e}".format(alpha_y, salpha_y)

    # This is broken lul
    # if not fromcache:
    #     write_cache(filename, data)

    return outstring, (opts["dataname"], data, paramsX, paramsY, (alpha_x, salpha_x), (alpha_y, salpha_y),opts["current"])

def process_position_data(filename, opts={}):
    outstring = ""

    data, fromcache = read_datfile(filename, opts["skiprows"])

    akbx, sakbx = get_alpha_kb_ratio(data, opts["T"], units=opts["convertX"], sunits=opts["uconvertX"])
    akby, sakby = get_alpha_kb_ratio(data, opts["T"], column=1, units=opts["convertY"], sunits=opts["uconvertY"])

    outstring += "\n--------------------\n"
    outstring += "Position dataset \"{}\"\n".format(opts["dataname"])
    outstring += "(\\alpha / k_B)_x = {0:1.3e} \\pm {1:1.3e}\n(\\alpha / k_B)_y = {2:1.3e} \\pm {3:1.3e}".format(akbx, sakbx, akby, sakby)

    if not fromcache:
        write_cache(filename, data)

    return outstring, (opts["dataname"], (akbx, sakbx), (akby, sakby))

##### EXECUTION #####
T = 293.15 * u.K           # Assumed to be room temperature

convertXdata1 = ((573.26e-3) * (u.V / u.micron))
convertYdata1 = ((550.22e-3) * (u.V / u.micron))
uconvXdata1   = (.186e-3 * (u.V / u.micron))
uconvYdata1   = (.223e-3 * (u.V / u.micron))
current1      = 200.0

convertXdata3 = ((683.89e-3) * (u.V / u.micron))
convertYdata3 = ((559.76e-3) * (u.V / u.micron))
uconvXdata3   = (.166e-3 * (u.V / u.micron))
uconvYdata3   = (.250e-3 * (u.V / u.micron))
current3      = 89.9

convertXdata4 = ((745.20e-3) * (u.V / u.micron))
convertYdata4 = ((703.16e-3) * (u.V / u.micron))
uconvXdata4   = (.15e-3 * (u.V / u.micron))
uconvYdata4   = (.68e-3 * (u.V / u.micron))
current4      = 200.0

convertXdata5 = ((1.09) * (u.V / u.micron))
convertYdata5 = ((1.00) * (u.V / u.micron))
uconvXdata5   = (1.0e-3 * (u.V / u.micron))
uconvYdata5   = (1.0e-3 * (u.V / u.micron))
current5      = 275.0

def main():
    freqdata = []
    posdata  = []
    kb  = []
    skb = []
    alphax = []
    alphay = []
    salphax = []
    salphay = []
    currents = []

    # By default, the pool executor only spins off 5 threads. This should be enough for us.
    with ProcessPoolExecutor(max_workers=5) as pool:
        futures = []

        futures.append(pool.submit(process_position_data, "../data/data1.dat", opts={
            "dataname": "Data 1",
            "skiprows": 2,
            "convertX": convertXdata1,
            "convertY": convertYdata1,
            "uconvertX": uconvXdata1,
            "uconvertY": uconvYdata1,
            "T": T
        }))

        futures.append(pool.submit(process_position_data, "../data/data2.dat", opts={
            "dataname": "Data 2",
            "skiprows": 2,
            "convertX": convertXdata1,
            "convertY": convertYdata1,
            "uconvertX": uconvXdata1,
            "uconvertY": uconvYdata1,
            "T": T
        }))

        futures.append(pool.submit(process_position_data, "../data/data3.dat", opts={
            "dataname": "Data 3",
            "skiprows": 2,
            "convertX": convertXdata3,
            "convertY": convertYdata3,
            "uconvertX": uconvXdata3,
            "uconvertY": uconvYdata3,
            "T": T
        }))

        futures.append(pool.submit(process_position_data, "../data/data4.dat", opts={
            "dataname": "Data 4",
            "skiprows": 2,
            "convertX": convertXdata4,
            "convertY": convertYdata4,
            "uconvertX": uconvXdata4,
            "uconvertY": uconvYdata4,
            "T": T
        }))

        futures.append(pool.submit(process_position_data, "../data/data5.dat", opts={
            "dataname": "Data 5",
            "skiprows": 2,
            "convertX": convertXdata5,
            "convertY": convertYdata5,
            "uconvertX": uconvXdata5,
            "uconvertY": uconvYdata5,
            "T": T
        }))

        futures.append(pool.submit(process_frequency_data, "../data/freq1.FDdat", opts={
            "dataname": "Data 1",
            "skiprows": 4,
            "T": T,
            "current": current1
        }))

        futures.append(pool.submit(process_frequency_data, "../data/freq3.FDdat", opts={
            "dataname": "Data 3",
            "skiprows": 4,
            "T": T,
            "current": current3
        }))

        futures.append(pool.submit(process_frequency_data, "../data/freq4.FDdat", opts={
            "dataname": "Data 4",
            "skiprows": 4,
            "T": T,
            "current": current4
        }))

        futures.append(pool.submit(process_frequency_data, "../data/freq5.FDdat", opts={
            "dataname": "Data 5",
            "skiprows": 4,
            "T": T,
            "current": current5
        }))

        for x in wait(futures)[0]:
            res = x.result()

            # Print out the output debugging statements from each result
            print(res[0])

            # Check if future is from freq space
            if len(res[1]) == 7:
                freqdata.append(res[1])
            else:
                posdata.append(res[1])

    for ftup in freqdata:
        # plot_power_spectrum(ftup)

        for ptup in posdata:
            if ptup[0] == ftup[0]:
                kbx, kby, skbx, skby, alpha_x, alpha_y, salpha_x, salpha_y= extract_kb_alpha(ptup, ftup)
                kb.append(kbx.value)
                skb.append(skbx.value)
                kb.append(kby.value)
                skb.append(skby.value)
                alphax.append(alpha_x.value)
                alphay.append(alpha_y.value)
                salphax.append(salpha_x.value)
                salphay.append(salpha_y.value)
                currents.append(ftup[6])
                #currents.append(ftup[6])

    kb  = np.array(kb)
    skb = np.array(skb)
    alphax = np.array(alphax)
    alphay = np.array(alphay)
    salphax = np.array(salphax)
    salphay = np.array(salphay)
    currents = np.array(currents)

    mu_kb = np.mean(kb)
    # s_kb  = np.sqrt(propagate(lambda a: np.mean(a), kb, skb)) <-- This is too narrow an uncertainty I think
    s_kb  = np.sum(skb)

    print("\n===== DETERMINED VALUE OF BOLTZMANN'S CONSTANT =====\n")
    print("kB = {0:1.3e} \\pm {1:1.3e} J / K\n".format(mu_kb, s_kb))
    plot_kb(mu_kb,s_kb)
    plot_current_dependence([alphax,alphay],[salphax,salphay],currents)
    plt.show()

if __name__ == "__main__":
    main()
