from dim4_pos import *


def energy(field_, kappa_):
    e = field_ * np.roll(field_, +1, 0)
    e += field_ * np.roll(field_, +1, 1)
    e += field_ * np.roll(field_, +1, 2)
    e += field_ * np.roll(field_, +1, 3)
    return -kappa_ * np.sum(e)


def order_param(field_):
    pp = np.sum((field_ + 1))/2
    return (pp - (v - pp))/v


def c2(field_, p_=(0, 0, 0)):
    px = np.cos(np.arange(nx) * p_[0])
    py = np.cos(np.arange(ny) * p_[1])
    pz = np.cos(np.arange(nz) * p_[2])
    tmp1 = np.sum(field_ * pz.reshape(1,nz,1,1), 1)
    tmp1 = np.sum(tmp1 * py.reshape(1,ny,1), 1)
    tmp1 = np.sum(tmp1 * px.reshape(1,nx), 1)
    px = np.cos(np.arange(nx) * -p_[0])
    py = np.cos(np.arange(ny) * -p_[1])
    pz = np.cos(np.arange(nz) * -p_[2])
    tmp2 = np.sum(field_ * pz.reshape(1,nz,1,1), 1)
    tmp2 = np.sum(tmp2 * py.reshape(1,ny,1), 1)
    tmp2 = np.sum(tmp2 * px.reshape(1,nx), 1)
    return tmp1[0] * tmp2
