import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)

m_sun = 1.99e30

data = np.loadtxt("out/data.csv", delimiter=",", skiprows=1)
r = data[:,0]
p = data[:,1]
m = data[:,2]
rho = data[:,3]
e = data[:,4]

fig = plt.figure()
plt.plot(r/1e3,p,color='black')
plt.grid()
plt.xlabel("radius $r$ (km)")
plt.ylabel("pressure $p(r)$")
plt.title("pressure against radius")
fig.tight_layout()
fig.savefig('out/pres.png')

fig = plt.figure()
plt.plot(r/1e3,m/m_sun,color='black')
plt.grid()
plt.xlabel("radius $r$ (km)")
plt.ylabel(r"mass $m(r)/m_\mathrm{sun}$")
plt.title("mass against radius")
fig.tight_layout()
fig.savefig('out/mass.png')

fig = plt.figure()
plt.plot(r/1e3,rho,color='black')
plt.grid()
plt.xlabel("radius $r$ (km)")
plt.ylabel(r"density $\rho(r)$")
plt.title("density against radius")
fig.tight_layout()
fig.savefig('out/dens.png')