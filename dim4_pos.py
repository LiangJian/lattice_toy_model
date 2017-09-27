import numpy as np

nt = 8
nz = 4
ny = 10
nx = 10

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
index = np.array(index)
