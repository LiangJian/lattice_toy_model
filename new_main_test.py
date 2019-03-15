from updates import *
from new_HK import *
from measure import energy, order_param, c2
import time

fields = start(2)  # two fields
phi = fields[0]
rho = fields[1]

kappa_eff = kappa_eff(rho)

st = time.time()
print('whole lattice', v)
c2s = []
nconf = 100
for i in range(nconf):
    print(i)
    update(phi, rho)
    print(energy(phi, kappa_phi), order_param(phi),'\t', energy(rho, kappa_phi), order_param(rho))
    c2s.append(c2(phi, p_=(1,0,0)))
ed = time.time()
c2s = np.array(c2s).reshape(nconf, nt)
c2s_ave = np.average(c2s, 0)
c2s_err = np.std(c2s, 0)/np.sqrt(nconf)
print(c2s_ave)
print(c2s_err)
print('done in %d s' % (ed - st))
print('==============================')
