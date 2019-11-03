import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)

m_sun = 1.99e30

data = np.loadtxt("out/mr.csv", delimiter=",", skiprows=1)
m = data[:,1]
r = data[:,2]

fig = plt.figure()
plt.plot(r/1e3,m/m_sun,color='black')
plt.grid()
plt.xlabel("radius $R$ (km)")
plt.ylabel("mass $M(R)/M_\mathrm{sun}$")
plt.title("$M$-$R$ curve of stable neutron stars")
fig.tight_layout()
fig.savefig('out/mr.png')
