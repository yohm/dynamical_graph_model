import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dat = np.loadtxt('lifetime.dat', delimiter=' ')

#taking linear binning
max_bin = 100
bin_size = math.ceil( dat[:,0].max() / max_bin )
d = {}
for i in range(max_bin):
    d[(i+1)*bin_size] = 0
for x,y in dat:
    if x > 0:
        k = bin_size*( math.ceil(float(x)/bin_size) )
        d[k] = d.get(k,0) + y
normalized = { k:float(v)/k for k,v in d.items() }
x = sorted( normalized.keys() )
y = [ normalized[k] for k in x ]

plt.xlabel( "lifetime" )
plt.ylabel( "frequency" )

plt.xscale('linear')
plt.yscale('log')
plt.plot( x, y, 'o-' )

plt.savefig("lifetime_linbin.png")
#plt.show()
