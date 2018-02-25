# pyTrans
Python fitter modified from pyExcal to transform between usno and the Kepler magnitudes

How to run pyTrans.

1.  Run the following command.
    ./pyTrans.py --inputFile keplerFinal.csv --band1 gU --band2 gM
    
2.  Usage.

usage: pyTrans2.py [-h] [--inputFile INPUTFILE] [--band1 BAND1]
                   [--band2 BAND2] [--verbose VERBOSE]

optional arguments:
  -h, --help            show this help message and exit
  --inputFile INPUTFILE
                        Name of input file
  --band1 BAND1         comma-separated list of filter bands to consider
  --band2 BAND2         comma-separated list of filter bands to consider
  --verbose VERBOSE     verbosity level of output to screen (0,1,2,...)

3.  You should get something similar to the following output.

 Let's Get This Pary Started.... 

inputFile =  keplerFinal.csv
Reference band =  gU
Band to be calculated =  gM

 Lets check on the status of things .... 
 
band1 =  gU   band2 =  gM   cband =  rM

mdn =  0.026

Initial parameter values:   [0.025999999999999801, 0.0]

Converged

Converged with chi squared  0.0525854897245 

degrees of freedom, dof  93 

RMS of residuals (i.e. sqrt(chisq/dof))  0.0237788850323

Reduced chisq (i.e. variance of residuals)  0.000565435373381

Fitted parameters at minimum, with 68% C.I.:

a               0.027529 +/-    0.00871078   (31.642225 percent)

k           -0.000111846 +/-     0.0122523   (10954.683792 percent)


Correlation matrix:
                a          k    
                
a            1.000000

k           -0.959978   1.000000


Your fit equation is:

   gU - gM = 0.028 + -0.00011*(gM-rM)

Outputting QA plot qa-keplerFinal.csv_color.gM-band.png


That's all, folks!


    
