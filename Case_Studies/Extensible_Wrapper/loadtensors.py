import numpy as np
import matplotlib.pyplot as plt
import scipy.io 

fschmidt = scipy.io.FortranFile ( "schmidttest.dat" )
schmidtdata = []
try:
    schmidtdata.append (fschmidt.read_reals()) 
except IOError:
    print 'done parsing schmidt values'
schmidtdata = np.array (schmidtdata)

gammadata = []
fgamma = scipy.io.FortranFile ( "gammatest.dat" )
try:
    gammadata.append (fgamma.read_reals()) 
except IOError:
    print 'done parsing gamma values'
gammadata = np.array (gammadata)

