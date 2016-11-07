
import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle


xlist = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
ylist = [1,2,3,4,5,6,7]

fig = plt.figure()
plt.plot( xlist, ylist, 'r-', markersize=4, markerfacecolor='r', linewidth=2., label="Test Data")

plt.gca().set_xticks([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])

plt.xlabel("Livetime [yr]")
plt.ylabel("T$_{1/2}$ [10$^{25}$ yr]")

plt.ylim([0, 10])
plt.legend(loc="lower right")
#plt.legend(loc="upper left")
fig.set_size_inches(6,5)

plt.savefig( "test_data.pdf" )

plt.show()

