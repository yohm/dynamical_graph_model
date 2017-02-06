import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dat = np.loadtxt("timeseries.dat", delimiter=' ')

plt.subplot(2,1,1)

plt.xlabel("Time (thousand steps)")
plt.ylabel("Number of species")
plt.plot( dat[:,0]/1000, dat[:,1], 'o-', label='# of species' )

plt.subplot(2,1,2)
plt.xlabel("Time (thousand steps)")
plt.ylabel("")
plt.plot( dat[:,0]/1000, dat[:,2], 'o-', label='link density' )
plt.plot( dat[:,0]/1000, dat[:,3], 'o-', label='CC' )
plt.legend(loc='best')

plt.tight_layout()
plt.savefig("timeseries.png")
#plt.show()
