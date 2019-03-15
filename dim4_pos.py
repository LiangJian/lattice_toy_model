import numpy as np

#  lattice size as global params
nt = 16
nz = 4
ny = 4
nx = 4

#  t z y x
ns = np.array((nt, nz, ny, nx))
vs = nx * ny * nz
v = vs * nt

sites = np.arange(0, v)

index = []
offset = np.zeros(shape=(nt, nz, ny, nx))
count = 0
for it in range(nt):
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                index.append((it, iz, iy, ix))
                offset[it, iz, iy, ix] = count


def move_back(index_, mu_):
    index_p_ = list(index_)
    index_p_[mu_] -= 1
    if (index_p_[mu_]) < 0:
        index_p_[mu_] += ns[mu_]
    return tuple(index_p_)
