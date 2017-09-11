##
# Hoshen–Kopelman algorithm for 4 dim
##


from dim4_pos import *
import numpy as np


lattice = np.random.randint(0, 2, v).reshape(nt, nz, ny, nx)
label = np.random.randint(0, 1, v).reshape(nt, nz, ny, nx)


def equive(pos_, label1_, label2_):
    for it_ in range(pos_[0] + 1):
        for iz_ in range(nz):
            for iy_ in range(ny):
                for ix_ in range(nx):
                    if label[it_, iz_, iy_, ix_] == label1_:
                        label[it_, iz_, iy_, ix_] = label2_


largest_label = 0
for ix in sites:
    index_ = index[ix]
    if lattice[tuple(index_)]:
        xm = np.zeros(4, dtype='int')
        for iu in range(4):
            if index_[iu] == 0:
                xm[iu] = 0
            else:
                index_[iu] -= 1
                xm[iu] = lattice[tuple(index_)]
                index_[iu] += 1

            # no bond
            if not xm[0] and not xm[1] and not xm[2] and not xm[3]:
                largest_label = largest_label + 1
                label[tuple(index_)] = largest_label

            # one bond
            # t
            elif xm[0] and not xm[1] and not xm[2] and not xm[3]:
                label[tuple(index_)] = label[index_[0]-1, index_[1], index_[2], index_[3]]
            # z
            elif not xm[0] and xm[1] and not xm[2] and not xm[3]:
                label[tuple(index_)] = label[index_[0], index_[1]-1, index_[2], index_[3]]
            # y
            elif not xm[0] and not xm[1] and xm[2] and not xm[3]:
                label[tuple(index_)] = label[index_[0], index_[1], index_[2]-1, index_[3]]
            # x
            elif not xm[0] and not xm[1] and not xm[2] and xm[3]:
                label[tuple(index_)] = label[index_[0], index_[1], index_[2], index_[3]-1]

            # two bonds
            # Since only the outermost loop may be scanned by part, order t z y x is somehow better
            # t z
            elif xm[0] and xm[1] and not xm[2] and not xm[3]:
                equive(tuple(index_),
                      label[index_[0]-1, index_[1], index_[2], index_[3]],
                      label[index_[0], index_[1]-1, index_[2], index_[3]])
                label[tuple(index_)] = label[index_[0], index_[1]-1, index_[2], index_[3]]

            # t y
            elif xm[0] and not xm[1] and xm[2] and not xm[3]:
                equive(tuple(index_),
                      label[index_[0]-1, index_[1], index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2]-1, index_[3]])
                label[tuple(index_)] = label[index_[0], index_[1], index_[2]-1, index_[3]]
            # t x
            elif xm[0] and not xm[1] and not xm[2] and xm[3]:
                equive(tuple(index_),
                      label[index_[0]-1, index_[1], index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]]-1)
                label[tuple(index_)] = label[index_[0], index_[1], index_[2], index_[3]]-1
            # z y
            elif not xm[0] and xm[1] and xm[2] and not xm[3]:
                equive(tuple(index_),
                      label[index_[0], index_[1]-1, index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2]-1, index_[3]])
                label[tuple(index_)] = label[index_[0], index_[1], index_[2]-1, index_[3]]
            # z x
            elif not xm[0] and xm[1] and not xm[2] and xm[3]:
                equive(tuple(index_),
                      label[index_[0], index_[1]-1, index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                label[tuple(index_)] = label[index_[0], index_[1], index_[2], index_[3]-1]
            # y x
            elif not xm[0] and not xm[1] and xm[2] and xm[3]:
                equive(tuple(index_),
                      label[index_[0], index_[1], index_[2]-1, index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                label[tuple(index_)] = label[index_[0], index_[1], index_[2], index_[3]-1]

            # three bonds
            # t z y
            elif xm[0] and xm[1] and xm[2] and not xm[3]:
                equive(tuple(index_),
                      label[index_[0]-1, index_[1], index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2]-1, index_[3]])
                equive(tuple(index_),
                      label[index_[0], index_[1]-1, index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2]-1, index_[3]])
                label[tuple(index_)] = label[index_[0], index_[1], index_[2]-1, index_[3]]
            # t z x
            elif xm[0] and xm[1] and not xm[2] and xm[3]:
                equive(tuple(index_),
                      label[index_[0]-1, index_[1], index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                equive(tuple(index_),
                      label[index_[0], index_[1]-1, index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                label[tuple(index_)] = label[index_[0], index_[1], index_[2], index_[3]-1]
            # t y x
            elif xm[0] and not xm[1] and xm[2] and xm[3]:
                equive(tuple(index_),
                      label[index_[0]-1, index_[1], index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                equive(tuple(index_),
                      label[index_[0], index_[1], index_[2]-1, index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                label[tuple(index_)] = label[index_[0], index_[1], index_[2], index_[3]-1]
            # z y x
            elif not xm[0] and xm[1] and xm[2] and xm[3]:
                equive(tuple(index_),
                      label[index_[0], index_[1]-1, index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                equive(tuple(index_),
                      label[index_[0], index_[1], index_[2]-1, index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                label[tuple(index_)] = label[index_[0], index_[1], index_[2], index_[3]-1]

            # four bonds
            # t z y x
            else:
                equive(tuple(index_),
                      label[index_[0]-1, index_[1], index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                equive(tuple(index_),
                      label[index_[0], index_[1]-1, index_[2], index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                equive(tuple(index_),
                      label[index_[0], index_[1], index_[2]-1, index_[3]],
                      label[index_[0], index_[1], index_[2], index_[3]-1])
                label[tuple(index_)] = label[index_[0], index_[1], index_[2], index_[3]-1]

print(lattice[..., 0, 0].reshape(nt, nz))
print(label[..., 0, 0].reshape(nt, nz))
