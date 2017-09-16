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


def prob_phi(kappa_eff_, phi_):
    tmp0_ = ((1 - np.exp(-2 * kappa_eff_[:, 0])) *
             ((phi_ / np.roll(phi_, -1, 0)) + 1) / 2).reshape((nt, ny, nz, nx, 1))
    tmp1_ = ((1 - np.exp(-2 * kappa_eff_[:, 1])) *
             ((phi_ / np.roll(phi_, -1, 1)) + 1) / 2).reshape((nt, ny, nz, nx, 1))
    tmp2_ = ((1 - np.exp(-2 * kappa_eff_[:, 2])) *
             ((phi_ / np.roll(phi_, -1, 2)) + 1) / 2).reshape((nt, ny, nz, nx, 1))
    tmp3_ = ((1 - np.exp(-2 * kappa_eff_[:, 3])) *
             ((phi_ / np.roll(phi_, -1, 3)) + 1) / 2).reshape((nt, ny, nz, nx, 1))
    return np.concatenate((tmp0_, tmp1_, tmp2_, tmp3_), 4)


def prob_rho(rho_, kappa_rho_):
    p_ = 1 - np.exp(-2 * kappa_rho_)
    tmp0_ = p_ * ((rho_ / np.roll(rho_, -1, 0) + 1) / 2).reshape(nt, ny, nz, nx, 1)
    tmp1_ = p_ * ((rho_ / np.roll(rho_, -1, 1) + 1) / 2).reshape(nt, ny, nz, nx, 1)
    tmp2_ = p_ * ((rho_ / np.roll(rho_, -1, 2) + 1) / 2).reshape(nt, ny, nz, nx, 1)
    tmp3_ = p_ * ((rho_ / np.roll(rho_, -1, 3) + 1) / 2).reshape(nt, ny, nz, nx, 1)
    return np.concatenate((tmp0_, tmp1_, tmp2_, tmp3_), 4)


def prob_rho_cluster_kernel(phi_, x_, mu_):
    pass


def prob_rho_cluster(phi_, x_, mu_):
    pass


def update(phi_, rho_):
    # step 1, using 'prob_phi' to generate phi field clusters
    #         Actually the process of generating clusters should be a modified HK process
    # step 2, after the clusters are formed, set the spin of the clusters with equal probabilities
    #         So it's better to have a index (list) of all clusters and also to have quick access to
    #         the points in each cluster (list of list)
    # step 3, using 'prob_rho' to generate rho field clusters
    # step 4, after the clusters are formed, set the spin of the clusters with 'prob_rho_cluster'
    pass