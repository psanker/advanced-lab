'''
This routine reads in and plots a sequence of spincore '.txt' files.  The plots 
are saved to disk with filename 'xxx.png', where 'xxx.txt' is the file name 
containing the data.  Derived from the program plot_ispin_forloop, written by 
Greg Lemberskiy.  

Last update:  1/30/2012, 10/7/2012, 
              1/21/2013 (changed .pdf to .png) by Tycho Sleator
'''

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

# The following is the filename containing the data
#fname = "sample_NMR_file.txt"
fnums = range(40)

for fname in fnums:    
    name  = "stim_echo_sweep_r5"+str(fname)+".txt"
    print name
    infile = open(name,"r")
    text = infile.read()      # read file into a string
    infile.close()

    index = text.find("@SW=") # Find the position of "@SW="
    text2 = text[index:-1]    # Create new string beginning at "@SW="
    index2 = text2.find('\n') # Find the next CR in this string

    bw = float(text2[4:index2]) # This gives the bandwidth
    print 'bw = ',bw
    print '1/bw = ',1/bw  # Note that the time interval between points is 1/bw
    
# Read the data from the the file starting on line 13
    s1 = mlab.csv2rec(name, skiprows=12)  
    
    t =  (1/bw)*np.arange(len(s1)/2)  #time data
    
# assign the data to variables with shorter names
    s = s1['data']
    rs = s.reshape(-1,2) 
    rtp= np.transpose(rs) # rtp[0] is the real part and rtp[1] is the 
                          # imaginary part of the data
    
    za = rtp[0] + rtp[1]*1j # create complex array
    
    
# get maximum value of the data

# create the figure
    fig1   = plt.figure(figsize=(8,10))
    ax1    = fig1.add_subplot(211)
    ax2    = fig1.add_subplot(212)


# Top Figure (Re and Im pargs): ax1
# draw x and y axes
    ax1.axhline(color ='k')
    ax1.axvline(color ='k')
    print 'len(az.real)=',len(za.real)
    print 'len(t)=',len(t)

# plot the points
    ax1.plot(t,za.real, '-b')  # plot the real part (blue)
    ax1.plot(t,za.imag, '-r')  # plot the imaginary part (red)

# label the axes
    ax1.set_xlabel('Time (sec)',fontsize=14)
    ax1.set_ylabel('Signal',fontsize=14)
# specify the plot limits
    ax1.set_xlim(t[0],t[-1])
    ax1.set_ylim(-1.2e8, 1.0e8)

# Bottom Figure (magnitude): ax2 
# draw x and y axes
    ax2.axhline(color ='k')
    ax2.axvline(color ='k')

    magnitude = ((za.real)**2 + (za.imag)**2)**(0.5)

# plot the points
    ax2.plot(t,magnitude)
# label the axes
    ax2.set_xlabel('Time (sec)',fontsize=14)
    ax2.set_ylabel('Signal',fontsize=14)
# specify the plot limits
    ax2.set_xlim(t[0],t[-1])
    ax2.set_ylim(0, 1.2e8)
# display the figure

    plt.savefig(name.replace(".txt",".png"))
    plt.show()
