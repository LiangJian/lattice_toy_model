#############################################################
# Swendsen-Wang type percolation cluster algorithms
# Hoshen–Kopelman algorithm for labeling clusters on a grid
#############################################################

from dim4_pos import *


def start(n_field_):
    fields = []
    for i in range(n_field_):
        fields.append(np.random.randint(0, 2, v).reshape(tuple(ns)))
        fields[i] = fields[i]*2 - 1
    return fields

kappa_phi = 0.07325
kappa_rho = 0.0718
g = 0.008


def kappa_eff(rho_):
    tmp0_ = kappa_phi - 0.5 * g * (rho_ + np.roll(rho_, -1, 0)).reshape((nt, ny, nz, nx, 1))
    tmp1_ = kappa_phi - 0.5 * g * (rho_ + np.roll(rho_, -1, 1)).reshape((nt, ny, nz, nx, 1))
    tmp2_ = kappa_phi - 0.5 * g * (rho_ + np.roll(rho_, -1, 2)).reshape((nt, ny, nz, nx, 1))
    tmp3_ = kappa_phi - 0.5 * g * (rho_ + np.roll(rho_, -1, 3)).reshape((nt, ny, nz, nx, 1))
    return np.concatenate((tmp0_, tmp1_, tmp2_, tmp3_), 4)


def prob_kappa():
    pass


def prob_rho():
    pass