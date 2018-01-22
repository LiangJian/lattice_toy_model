from updates import *
from measure import *
from Hk_labelling import *
import time

# spin configurations with random spin orientations
fields = start(2)
phi = fields[0]
rho = fields[1]
kappa_eff = kappa_eff(rho)

print('whole lattice', v)
labels = bond(phi, 1)
print(count(labels))
labels = bond(phi, -1)
print(count(labels))
print('==============================')

label_list, cluster_points = handle_labels(labels)
print(len(label_list))
for i in label_list:
    print(len(cluster_points[i]), end='\t')
print('')


print('==============================')
labels = bond_prop(phi, -1, np.zeros(shape=(nt, nz, ny, nx, 4)) + 0.5)
print(count(labels))
label_list, cluster_points = handle_labels(labels)
print(len(label_list))
for i in label_list:
    print(len(cluster_points[i]), end='\t')
print('')

st = time.time()
for i in range(50):
    update(phi, rho)
    # print(energy(phi, kappa_phi), energy(rho, kappa_rho))
    print(order_param(phi), order_param(rho))
ed = time.time()
print('done in %d s' % (ed - st))
