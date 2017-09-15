from updates import *
from Hk_labelling import *

# spin configurations with random spin orientations
fields = start(2)
phi = fields[0]
rho = fields[1]
kappa_eff = kappa_eff(rho)

labels = bond(phi, 1)
print(count(labels))
labels = bond(phi, -1)
print(count(labels))
print(v)