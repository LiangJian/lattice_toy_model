#############################################################
# Swendsen-Wang type percolation cluster algorithms
# Hoshen–Kopelman algorithm for labeling clusters on a grid
# https://arxiv.org/pdf/hep-lat/9503028.pdf
#############################################################

from dim4_pos import *
from HK_labelling_prop import hk4_2_prop, handle_labels
import copy
import numba


@numba.jit
def start(n_field_):
    fields = []
    for i in range(n_field_):
        fields.append(np.random.randint(0, 2, v).reshape(tuple(ns)))
        fields[i] = fields[i]*2 - 1  # plus 1 or minus 1
    return fields


kappa_phi = 0.04
kappa_rho = 0.04
g = 0.00


@numba.jit
def kappa_eff(rho_):
    tmp0_ = kappa_phi - 0.5 * g * (rho_ + np.roll(rho_, -1, 0)).reshape((nt, nz, ny, nx, 1))
    tmp1_ = kappa_phi - 0.5 * g * (rho_ + np.roll(rho_, -1, 1)).reshape((nt, nz, ny, nx, 1))
    tmp2_ = kappa_phi - 0.5 * g * (rho_ + np.roll(rho_, -1, 2)).reshape((nt, nz, ny, nx, 1))
    tmp3_ = kappa_phi - 0.5 * g * (rho_ + np.roll(rho_, -1, 3)).reshape((nt, nz, ny, nx, 1))
    return np.concatenate((tmp0_, tmp1_, tmp2_, tmp3_), 4)


@numba.jit
def prob_phi(kappa_eff_):
    dice = np.random.random((nt, nz, ny, nx, 4))
    return dice < 1 - np.exp(-2 * kappa_eff_)

# dice = np.random.random((nt, nz, ny, nx, 4))

@numba.jit
def prob_rho(kappa_rho_):
    dice = np.random.random((nt, nz, ny, nx, 4))
    return dice < 1 - np.exp(-2 * kappa_rho_)


@numba.jit
def update(phi_, rho_):

    # step 1, using 'prob_phi' to generate phi field clusters
    #         Actually the process of generating clusters should be a modified HK process
    # step 2, after the clusters are formed, set the spin of the clusters with equal probabilities
    #         So it's better to have a index (list) of all clusters and also to have quick access to
    #         the points in each cluster (list of list)
    # step 3, using 'prob_rho' to generate rho field clusters
    # step 4, after the clusters are formed, set the spin of the clusters with 'prob_rho_cluster'

    labels_1_ = phi_.copy() * 0.0
    labels_2_ = phi_.copy() * 0.0
    hk4_2_prop(phi_, labels_1_, +1, prob_phi(kappa_eff(rho_)))
    hk4_2_prop(phi_, labels_2_, -1, prob_phi(kappa_eff(rho_)))

    label_list, cluster_points = handle_labels(labels_1_)
    for i in label_list:
        #print('len', len(cluster_points[i]))
        if np.random.random() < 0.5:
            for j in range(len(cluster_points[i])):
                phi_[cluster_points[i][j]] = -1

    label_list, cluster_points = handle_labels(labels_2_)
    for i in label_list:
        if np.random.random() < 0.5:
            for j in range(len(cluster_points[i])):
                phi_[cluster_points[i][j]] = +1


    labels_1_ *= 0.0
    labels_2_ *= 0.0
    hk4_2_prop(rho_, labels_1_, +1, prob_rho(kappa_rho))
    hk4_2_prop(rho_, labels_2_, -1, prob_rho(kappa_rho))

    label_list, cluster_points = handle_labels(labels_1_)
    for i in label_list:
        summation = 0
        for j in range(len(cluster_points[i])):
            for k in range(4):
                index0 = list(cluster_points[i][j])
                indexm = copy.deepcopy(index0)
                indexm[k] = (indexm[k] - 1) % ns[k]
                indexp = copy.deepcopy(index0)
                indexp[k] = (indexp[k] + 1) % ns[k]
                summation += phi_[tuple(index0)] * (phi_[tuple(indexm)] + phi_[tuple(indexp)])
        prop_cluster = np.exp(g * summation)
        prop_cluster = 1. / (1 + prop_cluster)
        if np.random.random() >= prop_cluster:
            for j in range(len(cluster_points[i])):
                rho_[cluster_points[i][j]] = -1

    label_list, cluster_points = handle_labels(labels_2_)
    for i in label_list:
        summation = 0
        for j in range(len(cluster_points[i])):
            for k in range(4):
                index0 = list(cluster_points[i][j])
                indexm = copy.deepcopy(index0)
                indexm[k] = (indexm[k] - 1) % ns[k]
                indexp = copy.deepcopy(index0)
                indexp[k] = (indexp[k] + 1) % ns[k]
                summation += phi_[tuple(index0)] * (phi_[tuple(indexm)] + phi_[tuple(indexp)])
        prop_cluster = np.exp(g * summation)
        prop_cluster = 1. / (1 + prop_cluster)
        if np.random.random() < prop_cluster:
            for j in range(len(cluster_points[i])):
                rho_[cluster_points[i][j]] = +1

