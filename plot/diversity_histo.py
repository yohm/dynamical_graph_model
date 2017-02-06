import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dat = np.loadtxt('diversity_histo.dat', delimiter=' ')

plt.xlabel( "number of species" )
plt.ylabel( "frequency" )

plt.yscale('log')
plt.plot( dat[:,0], dat[:,1], 'o-' )

plt.savefig("diversity_histo.png")
#plt.show()
