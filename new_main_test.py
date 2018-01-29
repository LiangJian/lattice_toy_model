from updates import *
from new_HK import *
import time

fields = start(2)
phi = fields[0]
rho = fields[1]
kappa_eff = kappa_eff(rho)

st = time.time()
print('whole lattice', v)
labels, labels_list = bond_prop_new(phi, 1, prob_phi(kappa_eff, phi))
print(labels)
ed = time.time()
print('done in %d s' % (ed - st))
print('==============================')