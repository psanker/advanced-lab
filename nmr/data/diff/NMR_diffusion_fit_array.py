'''
This routine to uses the Levenburg-Marquardt algorithm to fit a set of data
points to the function A*np.exp(-((gamma*t)**3)/6) + y0, which describes the
amplitude of a Hahn echo in the presence of a magnetic field gradient and
diffusion.  From the fit, one can determine the constant of self diffusion.
The Hahn echo amplitudes are obtained from an array pasted into the program.
Formally called NMR_diffusion_fit_sec

Last Update 9/29/2012, 10/7/2012 at 12:31 pm by Tycho Sleator
Amended by Patrick Anker & Kaitlyn Morrell (02/28/2018)
'''
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
# PK: Removed old imports

# Read the data from the .csv file "data.csv"
#s1 = mlab.csv2rec("Diffusion_data.csv", skiprows=2)

# print out the data names (and type) read in from csv file
#print(s1.dtype)

# assign the data to variables with shorter names
#t = s1.time_ms
#s = s1.echo_size
#sigs = s1.error
print("\n\n")
print("---------------------------------------------------------------------------")
print("File NMR_diffusion_fit_2.py")
print("This routine to uses the Levenburg-Marquardt algorithm to fit a set of data")
print("points to the function A*np.exp(-((gamma*t)**3)/6) + y0, which describes the")
print("amplitude of a Hahn echo in the presence of a magnetic field gradient and")
print("diffusion.  From the fit, one can determine the constant of self diffusion.")
print("---------------------------------------------------------------------------")
# PK: Converted python2 syntax to futureproof python3 syntax

t = 0.001 * np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140])
s =  np.array([0.74373218544388053, 0.54829865107798392, 0.82055123499403715, 1.0, 0.94426665226231288, 0.92975789160515865, 0.90893941276699119, 0.92556506362147972, 0.91820947988573054, 0.85451130141362652, 0.8240415919133246, 0.81537443267987475, 0.81178126883721502, 0.7646432792463588])
gamma = 5.59 #gyromagnetic ratio
print("t = ", t)
print("s = ", s)


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
    return A*np.exp(-((gamma**2)*alpha*(t**3))/12) + y0

# Initial guesses: modify these if fit doesn't converge.
A0 = 1.0
alpha0 = 10**(3)
y00 = 0.0

# This is the function that does the nonlinear fit:
#nlfit, nlpcov = scipy.optimize.curve_fit(f, t, s, p0=[A0, gamma0, y00], sigma=None)
nlfit, nlpcov = scipy.optimize.curve_fit(f, t, s, p0=[A0, alpha0, y00], sigma=None)

# These are the parameter estimates (assigned to more readible names):
A_nlfit = nlfit[0]
alpha_nlfit = nlfit[1]
y0_nlfit = nlfit[2]
print('A = {}'.format( A_nlfit))
print('alpha = {}'.format( alpha_nlfit))
print('y0 = {}'.format( y0_nlfit))

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
plt.figure(1, figsize = (8,4.5) )

# Draw x and y axes:
plt.axhline(color ='k')
plt.axvline(color ='k')

# Plot the points:
l1 = plt.plot(t, s, 'ob', label='Data')
# Include error bars:
#plt.errorbar(t, s, sigs, fmt='bo')

# Plot the fit:
l2 = plt.plot(tfit,sfit,'r', label='Fit')

# Plot the legend:
plt.legend(loc='upper right' )

# Plot some text with model equation and values of fit parameters.
# Use Latex equation style - $ ... $:
plt.text(0.010, 0.4,
    '$y \, =\, A e^{-(\\gamma^2\\alpha t^3)/12} + y_0$',
    ha='left', va='bottom', size='large')
plt.text(0.010, 0.3,
    '$\\alpha \, =\, ({0:6.4g} \pm {1:5.2g})/s^3$'.format(alpha_nlfit, alphasig_nlfit),
    ha='left', va='bottom', size='large')
plt.text(0.010, 0.2,
    '$A  \, =\, {0:6.3g} \pm  {1:6.1g}$'.format(A_nlfit, Asig_nlfit),
    ha='left', va='bottom', size='large')
plt.text(0.010, 0.1,
    '$y_0 \, =\, {0:6.3g} \pm {1:6.1g}$'.format(y0_nlfit, y0sig_nlfit),
    ha='left', va='bottom', size='large')


# Label the axes:
plt.xlabel('Time (s)',fontsize=14)
plt.ylabel('Echo Amplitude',fontsize=14)

# Specify the plot limits:
#plt.xlim(0.0,max(t))
plt.xlim(0.0,t[-1] + (t[1]-t[0])/2)
plt.ylim(-0.05,1.05)

# Display the figure:
plt.show()
