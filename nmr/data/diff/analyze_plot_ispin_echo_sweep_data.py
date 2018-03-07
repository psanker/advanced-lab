'''
File:  'analyze_plot_ispin_echo_sweep_data.py'
This file is derived from the file 'manipulate_plot_ispin_nmr_data.py'.
(which was previously updated:  2/20/2012 by Tycho Sleator)

This file reads in a sequence of spincore ".txt" files generated by the
"hahn_echo" sequence, extracts numbers from each of those files, generates
and prints an array from these numbers and plots the result.

The array can be pasted into a fitting program, such as
'NMR_diffusion_fit_array.py' for further analysis.

Last update:  10/07/2012 at 12:47 pm by Tycho Sleator
Amended by Patrick Anker & Kaitlyn Morrell (02/28/2018)
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.integrate as integrate
import scipy.optimize as opt

# PK: Changed file handling for more flexibility
import sys
import os
import glob
# PK: Removed unused matplotlib and scipy imports

'''
THE FOLLOWING ARE FUNCTION DEFINITIONS
'''
'''
Array y is considered a function of array x
Integrate y as function of x in the range of x1 to x2,
'''
def intsimps(y, x, x1, x2):
    x0 = x[0]

    dx = x[1] - x[0]
    # PK: Ensure datatype is int for indexing; np.floor() does not return int
    i1 = int(np.floor((x1-x0)/dx))  # index of starting value of x
    i2 = int(np.floor((x2-x0)/dx))  # index of ending value of x

    return integrate.simps(y[i1:i2], None, dx)


"""
This reads in complex data from a file and outputs a tuple (t,za),
where t is the time data and za is an array of complex numbers
corresponding to the time data
"""
def read_data_file(fname):
    infile = open(fname, "r")
    text = infile.read()      # read file into a string
    infile.close()

    index = text.find("@SW=") # Find the position of "@SW="
    text2 = text[index:-1]    # Create new string beginning at "@SW="
    index2 = text2.find('\n') # Find the next CR in this string
    # print 'text=',text2[0:index2]
    bw = float(text2[4:index2]) # This gives the bandwidth
    # print 'bw = ',bw
    # print '1/bw = ',1/bw  # Note that the time interval between points is 1/bw

    # Read the data from the the file starting on line 13
    s1 = mlab.csv2rec(fname, skiprows=12)
    npts = len(s1)/2  # number of complex data points

    # print 'npts = ',npts

    t =  (1/bw)*np.arange(npts) # time data

    # assign the data to variables with shorter names
    s = s1['data']
    rs = s.reshape(-1, 2)
    rtp = np.transpose(rs) # rtp[0] is the real part and rtp[1] is the imaginary
                           # part of the data

    return (t, (rtp[0] + rtp[1]*1j)) # create complex array

# rotate data in complex plane
#theta = -0.0*np.pi
# za = np.exp(1j*theta)*(rtp[0] + rtp[1]*1j) # create complex array

'''
EXECUTION BEGINS HERE
'''
def help():
    return "Usage: python <script name> [fit] <directory name>\nIf you want to fit the data from a directory, include the word \"fit\" in one of the arguments after the script file name"

# PK: Changed from file-based search to directory-based search for different datasets
def main(argv):
    fit  = False
    dirs = [] # Directories to process

    if len(argv) == 0:
        print("No directory given")
        print(help())
        return
    else:
        for arg in argv:
            path = os.path.join(os.getcwd(), arg)

            # Check if "fit" is provided
            if arg == "fit":
                print("Will fit data...")
                fit = True

            # Ensure directory exists
            elif not os.path.exists(path):
                print("Directory \"{}\" does not exist.".format(arg))
                continue
            else:
                dirs.append(path)

        for _dir in dirs:
            process_directory(sorted(glob.glob("{}/*.txt".format(_dir))), fit, _dir)

        plt.show() # Draw all plots

def process_directory(files, fit, dirname):
    taumin    = 10  # Minimum value of tau (ms)
    N_taus    = 14  # Number of taus to process
    Delta_tau = 10  # Step increase of tau (ms)

    taus  = taumin + Delta_tau * np.arange(N_taus)
    files = files[:N_taus]  # Filter the files to the number required

    echos = []

    for fname in files:
        # read data from file
        print("Reading: {}".format(fname))
        (t, za) = read_data_file(fname)

        echosize = intsimps(abs(za),t,0.017,0.023) # compute the 'size' of the echo
        echos.append(echosize)  # add to the array of echo sizes
        print('echosize =', echosize)

    rechos = echos/(max(echos))
    print('taus = ', list(taus))
    print('echo sizes = ', list(echos))
    print('relative echo sizes = ', list(rechos))

    '''
    CREATE THE FIGURE
    '''
    plt.figure(figsize=(8, 5))
    # draw x and y axes
    plt.axhline(color='k')
    plt.axvline(color='k')
    #
    # plot the points
    plt.plot(taus, rechos, 'ob')  # plot the real part (blue)
    #
    ## label the axes
    plt.xlabel('Echo Delay $\\tau$ (ms)', fontsize=14)
    plt.ylabel('Relative Echo Size', fontsize=14)
    #
    # specify the plot limits
    plt.ylim(0, 1.05)
    #

    if fit:
        print("Fitting dir \"{}\"...".format(dirname))
        fit_data(taus, rechos, dirname)

def fit_data(taus, rels, fname):
    t = 1e-3 * np.array(taus)
    s = rels

    '''
    -- information on nonlinear curve fitting --

    nlfit,nlpcov = scipy.optimize.curve_fit(f, xdata, ydata, p0=[a0, b0, ...],
      sigma=None, **kwargs)

    * f(xdata, a, b, ...): is the fitting function where xdata is an array of values
       of the independent variable and a, b, ... are the fitting parameters, however
       many there are, listed as separate arguments.  f(xdata, a, b, ...) should
       return the y value of the fitting function.
    * xdata: is the array containing the x data.
    * ydata: is the array containing the y data.
    * p0: is a tuple containing the initial guesses for the fitting parameters. The
      guesses for the fitting parameters are set equal to 1 if they are left
      unspecified.
    * sigma: is the array containing the uncertainties in the y data.
    * **kwargs: are keyword arguments that can be passed to the fitting routine
      numpy.optimize.leastsq that curve_fit calls. These are usually left unspecified.

    '''
    # Define nonlinear fitting function.
    def f(t, A, alpha, y0):
        return A*np.exp(-alpha*(t**3)) + y0

    # Initial guesses: modify these if fit doesn't converge.
    A0 = 1.0
    alpha0 = 10**(3)
    y00 = 0.0

    # This is the function that does the nonlinear fit:
    #nlfit, nlpcov = scipy.optimize.curve_fit(f, t, s, p0=[A0, gamma0, y00], sigma=None)
    nlfit, nlpcov = opt.curve_fit(f, t, s, p0=[A0, alpha0, y00], sigma=None)

    # These are the parameter estimates (assigned to more readible names):
    A_nlfit = nlfit[0]
    alpha_nlfit = nlfit[1]
    y0_nlfit = nlfit[2]

    '''
    Below are the uncertainties in the estimates:  (note that "nlpcov" is the
    "covariance matrix".  the diagonal elements of the covariance matrix (nlpcov[i][i])
    are the variances of the fit parameters, whose square roots give the standard
    deviation or uncertainties in the fit parameters.  The off-diagonal elements give
    the correlations between fit parameters.
    '''
    Asig_nlfit =  np.sqrt(nlpcov[0][0])
    alphasig_nlfit = np.sqrt(nlpcov[1][1])
    y0sig_nlfit = np.sqrt(nlpcov[2][2])

    '''
    The following statements generate formatted output of the fit parameters.
    The number before the colon refers to the sequence of arguments in the .format
    command.  In "{0:6.3f}", the 0 refers to the first argument in the .format()
    command, the 6 refers to the total width of the field, and the 3 refers to the
    number of digits after the decimal point.
    See http://docs.python.org/tutorial/inputoutput.html for more information.
    '''

    # Compute reduced chi square for nonlinear fit.
    #redchisqrNL = (((s - f(t, A_nlfit, alpha_nlfit, y0_nlfit))/sigs)**2).sum()/float(len(t)-2)

    # Create an array of times t for the purpose of plotting the fit:
    tfit = np.arange(0.0, t[-1] + (t[1]-t[0])/2, 0.0001)

    # Generate y-values from parameters determined from the curve fit, so that
    # we can plot the fit on top of the data:
    sfit = f(tfit, A_nlfit, alpha_nlfit, y0_nlfit)

    # Create the figure:
    fig = plt.figure(1, figsize = (8,4.5) )

    # Draw x and y axes:
    plt.axhline(color ='k')
    plt.axvline(color ='k')

    # Plot the points:
    plt.plot(t, s, 'ob', label='Data')
    # Include error bars:
    #plt.errorbar(t, s, sigs, fmt='bo')

    # Plot the fit:
    plt.plot(tfit,sfit,'r', label='Fit')

    # Plot the legend:
    plt.legend(loc='upper right' )

    # Plot some text with model equation and values of fit parameters.
    # Use Latex equation style - $ ... $:
    plt.text(0.010, 0.4, '$y \, =\, A e^{-\\alpha t^3} + y_0$', ha='left', va='bottom', size='large')
    plt.text(0.010, 0.3, '$\\alpha \, =\, ({0:6.4g} \pm {1:5.2g})/s^3$'.format(alpha_nlfit, alphasig_nlfit), ha='left', va='bottom', size='large')
    plt.text(0.010, 0.2, '$A  \, =\, {0:6.3g} \pm  {1:6.1g}$'.format(A_nlfit, Asig_nlfit), ha='left', va='bottom', size='large')
    plt.text(0.010, 0.1, '$y_0 \, =\, {0:6.3g} \pm {1:6.1g}$'.format(y0_nlfit, y0sig_nlfit), ha='left', va='bottom', size='large')


    # Label the axes:
    plt.xlabel('Time (s)', fontsize=14)
    plt.ylabel('Echo Amplitude', fontsize=14)

    # Specify the plot limits:
    #plt.xlim(0.0,max(t))
    plt.xlim(0.0,t[-1] + (t[1]-t[0])/2)
    plt.ylim(-0.05,1.05)

    # Save figure
    plt.savefig(fname)
    plt.close(fig)

if __name__ == "__main__":
    main(sys.argv[1:])
