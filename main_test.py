from updates import *
from Hk_labelling import *

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
labels = bond_prop(phi, -1, 0.5)
print(count(labels))
label_list, cluster_points = handle_labels(labels)
print(len(label_list))
for i in label_list:
    print(len(cluster_points[i]), end='\t')
print('')
