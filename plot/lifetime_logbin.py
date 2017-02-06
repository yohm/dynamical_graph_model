import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dat = np.loadtxt('lifetime.dat', delimiter=' ')

#taking log binning
d = {}
for x,y in dat:
    if x > 0:
        k = 2**( math.ceil(math.log2(x)) )
        d[k] = d.get(k,0) + y
normalized = { k:float(v)/k for k,v in d.items() }
x = sorted( normalized.keys() )
y = [ normalized[k] for k in x ]

plt.xlabel( "lifetime" )
plt.ylabel( "frequency" )

plt.xscale('log')
plt.yscale('log')
plt.plot( x, y, 'o-' )

plt.savefig("lifetime_logbin.png")
#plt.show()
