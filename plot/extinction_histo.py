import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dat = np.loadtxt('extinction_histo.dat', delimiter=' ')

plt.xlabel( "extinction size" )
plt.ylabel( "frequency" )

plt.yscale('log')
plt.plot( dat[:,0], dat[:,1], 'o-' )

plt.savefig("extinction_histo.png")
