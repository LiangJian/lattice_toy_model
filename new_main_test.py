from HK_labelling_prop import hk4_prop, hk4_2_prop
from updates import start
import time
from dim4_pos import *
from updates import update, kappa_phi, kappa_rho
from measure import c2, order_param, energy

fields = start(2)  # two fields
phi = fields[0]
rho = fields[1]
rho = phi.copy()

label = phi.copy()
props = np.zeros(shape=(nt, ny, nz, nx, 4)) == 0

st = time.time()
hk4_prop(phi, label, 1, props)
ed = time.time()
print('%.2f' % (ed - st))
bak = label.copy()

st = time.time()
hk4_2_prop(phi, label, 1, props)
ed = time.time()
print('%.2f' % (ed - st))
print(np.sum(bak - label))

# Tc of 4D case, 6.682(2)
# k = 1/6.682 = 0.1497
# 0.1497/2 ~ 0.07485

phi_p = []
rho_p = []
for i in range(100):
    print(i)
    update(rho, phi)
    #print('%.2f'%energy(rho, kappa_rho), '%.2f'%energy(phi, kappa_phi))
    #print('%.2f' % order_param(rho), '%.2f' % order_param(phi))
    if i > 9:
        phi_p.append(np.sum(phi, axis=(1, 2, 3)))
        rho_p.append(np.sum(rho, axis=(1, 2, 3)))

phi_p = np.array(phi_p)
rho_p = np.array(rho_p)
print(phi_p.shape)

phi_p_a = np.average(phi_p, 0)
phi_p = phi_p - phi_p_a

rho_p_a = np.average(rho_p, 0)
rho_p = rho_p - rho_p_a

print(np.average(phi_p[0]*phi_p, 0))
print(np.average(rho_p[0]*rho_p, 0))