from dim4_pos import *


def energy(field_, kappa_):
    E = field_ * np.roll(field_, +1, 0)
    E += field_ * np.roll(field_, +1, 1)
    E += field_ * np.roll(field_, +1, 2)
    E += field_ * np.roll(field_, +1, 3)
    return -kappa_ * np.sum(E)


def order_param(field_):
    pp = np.sum((field_ + 1))/2
    return (pp - (v - pp))/v
