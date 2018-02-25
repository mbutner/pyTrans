#!/usr/bin/env python

# pyTrans.py 
# A very basic fitter to the SDSS primed system photometric equations of the form: 
# 	u_usno - u_mddho = C0 + C1 * (u_mho - g_mho) 
# 	g_usno - g_mho = C0 + C1 * (g_mho - r_mho) 
# 	r_usno - r_mho = C0 + C1 * (g_mho - r_mho) 
# 	i_usno - i_mho = C0 + C1 * (r_mho - i_mho) 
# 	z_usno - z_mho = C0 + C1 * (i_mho - z_mho) 
#	
#	 Reads in a CSV file with the following columns and fits the data to
#    the relevant photometric equation (u',g',r',i',z'): 
#      Name,ra,dec, uU, gU, rU, iU, zU, uUe, gUe, rUe, iUe, zUe, uM, gM, rM, iM, zM  
#    (order is not important, but spelling and case are important)   
#    Example:
#    
#    pyTrans.py --help
#    pyTrans.py --inputFile std-rguiz-test.g.csv --band g --verbose 1
#
#####################################################################################3
print   
print " Let's Get This Pary Started.... " 
print
##################################


def main():

    import argparse
        
    supportedBandList = ['gU','rU','iU','zU','gM','rM','iM','zM']

    
# Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputFile', help='Name of input file', default='keplerFinal.csv')
    parser.add_argument('--band1', help='comma-separated list of filter bands to consider', default='gU,rU,iU,zU,gM,rM,iM,zM')
    parser.add_argument('--band2', help='comma-separated list of filter bands to consider', default='gU,rU,iU,zU,gM,rM,iM,zM')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    if args.band1 not in supportedBandList:
        print """Filter band %s is not found in the list of supported bands...""" % (band1)
        print """The list of supported filter bands is %s""" % (supportedBandList)
        print """(Note that the filter band names in this list are case-sensitive.)"""
        print """Exiting now!..."""
        print 
        return 1
    if args.band2 not in supportedBandList:
        print """Filter band %s is not found in the list of supported bands...""" % (band2)
        print """The list of supported filter bands is %s""" % (supportedBandList)
        print """(Note that the filter band names in this list are case-sensitive.)"""
        print """Exiting now!..."""
        print 
        return 1
 
    status = pyTrans_fit(args)
 
    print
    print
    print "That's all, folks!"
    print

    return 0


##################################

def pyTrans_fit(args):

    # Based on a scripts at
    # http://linuxgazette.net/115/andreasen.html (by Anders Andreasen)
    # and at
    # http://www.phy.uct.ac.za/courses/python/examples/fitresonance.py (University of Cape Town)

    import numpy as np
    import math
    import os
    import sys
    from scipy.optimize import leastsq
    import matplotlib.pyplot as plt

    if args.verbose>0:
        print
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'pyTrans_fit'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print

    inputFile = args.inputFile
    band1 = args.band1
    band2 = args.band2

    # Check for the existence of the input file...
    if os.path.isfile(inputFile)==False:
        print """File %s does not exist...""" % (inputFile)
        print """Exiting now!"""
        print
        return 1

     # Verify input data is correct
    print "inputFile = " , inputFile
    print "Reference band = ", band1
    print "Band to be calculated = ", band2  

    # Extract file basename (to be used to name output qa files...)
    baseName = os.path.basename(inputFile)
    
    # Read whole file into numpy arrays using np.genfromtxt...
    data = np.genfromtxt(inputFile,dtype=None,delimiter=',',names=True)

    # Python dictionary of conjugate filter bands...
    cbandDict = {'gU':'rU','rU':'iU','iU':'rU','zU':'iU','gM':'rM','rM':'iM','iM':'rM','zM':'iM'}

    # Python dictionary of color names...
    colorNameDict = {'gU':'gU-rU','rU':'rU-iU','iU':'rU-iU','zU':'iU-zU', 'gM':'gM-rM','rM':'rM-iM','iM':'rM-iM','zM':'iM-zM'}
    

    # Grab the correct conjugate filter band from the
    #  conjugate filter band Python dictionary...
    cband = cbandDict[band2]

    # Extract relevant numpy arrays for fit....
    dmag = data[band2] - data[band1]
    color = data[band2]-data[cband]

    print 
    print " Lets check on the status of things .... "
    print "band1 = ", band1 , " ", "band2 = ", band2, " ", "cband = ", cband
#   print "dmag = ", dmag 
#   print "color = ", color
    print

    # Calculate the median of dmag for use as an initial guess
    # for the overall zeropoint offset..
    mdn = np.median( dmag, None )
    print "mdn = ", mdn
    
# Parameter names
    pname = (['a', 'k'])

    # Initial parameter values
    p0 = [mdn, 0.0]

    print
    print 'Initial parameter values:  ', p0
     
   #dmag = str(dmag)
   #color = str(color)

   # Perform fit
    p,cov,infodict,mesg,ier = leastsq(residuals, p0, args=(color, dmag), maxfev=10000, full_output=1)
    if (ier >=1 and ier <=4):
        print "Converged"
    else:
        print "Not converged"
        print mesg

    # Calculate some descriptors of the fit 
    # (similar to the output from gnuplot 2d fits)
    chisq=sum(infodict['fvec']*infodict['fvec'])
    dof=len(dmag)-len(p)
    print "Converged with chi squared ",chisq , type(chisq)
    print "degrees of freedom, dof ", dof, type(dof)
    print "RMS of residuals (i.e. sqrt(chisq/dof)) ", math.sqrt(chisq/dof)
    print "Reduced chisq (i.e. variance of residuals) ", chisq/dof
    print

    # uncertainties are calculated as per gnuplot, "fixing" the result
    # for non unit values of the reduced chisq.
    # values at min match gnuplot
    print "Fitted parameters at minimum, with 68% C.I.:"
    for i,pmin in enumerate(p):
	print "%-10s %13g +/- %13g   (%5f percent)" % (pname[i],pmin,math.sqrt(cov[i,i])*math.sqrt(chisq/dof),100.*math.sqrt(cov[i,i])*math.sqrt(chisq/dof)/abs(pmin))
    print  

    print
    print "Correlation matrix:"
    # correlation matrix close to gnuplot
    print "               ",
    for i in range(len(pname)): print "%-10s" % (pname[i],),
    print
    for i in range(len(p)):
        print "%-10s" % pname[i],
        for j in range(i+1):
            print "%10f" % (cov[i,j]/math.sqrt(cov[i,i]*cov[j,j]),),
        #endfor
        print
    #endfor

    title="""%s - %s = %.3f + %.5f*(%s)""" % (band1, band2, p[0], p[1], colorNameDict[band2] )
    xlabel = ''
    ylabel = """%s - %s""" % (band2, band1)

    print
    print
    print """Your fit equation is:\n   %s""" % (title)
    print  
    # output QA plot...
#   qaPlot1 = """qa-%s_airmass.%s-band.png""" % (baseName, band)
#   print """Outputting QA plot %s""" % (qaPlot1)
#   xlabel = 'airmass'
#   plt.title(title)
#   plt.xlabel(xlabel)
#   plt.ylabel(ylabel)
#   plt.scatter(X,dmag)
#   colorMean=np.zeros(color.size)+np.mean(color)
#   plt.plot(X, fp(p,colorMean), '-', linewidth=2)
#   plt.grid(True)
#   plt.savefig(qaPlot1)
#   plt.clf()
  
    qaPlot2 = """qa-%s_color.%s-band.png""" % (baseName, band2)
    print """Outputting QA plot %s""" % (qaPlot2)
    xlabel = """%s - %s""" % (colorNameDict[band2], colorNameDict[band1])
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.scatter(color,dmag)
#   XMean=np.zeros(X.size)+np.mean(X)
    plt.plot(color, fp(p,color), 'r-', linewidth=2)
    plt.grid(True)
    plt.savefig(qaPlot2)
    plt.clf()

    return 0

##################################

# Parametric function:
#  p is the parameter vector;
#  color is the stars mho color

def fp(p,color):
    return p[0] + p[1]*color

##################################

# Error function:
def residuals(p,color,dmag):
    err = (dmag-fp(p,color))
    return err

##################################

if __name__ == "__main__":
    main()

##################################




