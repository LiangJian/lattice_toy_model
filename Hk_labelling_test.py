import numpy as np

'''
To test Hoshen–Kopelman algorithm
'''

nx = 8
ny = 8
nz = 8
nt = 8


def union2(label_, y_, label1_, label2_):
    for iy in range(y_ + 1):
        for ix in range(nx):
            if label_[iy, ix] == label1_:
                label_[iy, ix] = label2_


def hk2(lattice_, label_):
    largest_label = 0
    for y in range(ny):
        for x in range(nx):
            if lattice_[y, x] == 1:
                # Can I use other BD condition?
                if x == 0:
                    xm = 0
                else:
                    xm = lattice_[y, x-1]
                if y == 0:
                    ym = 0
                else:
                    ym = lattice_[y-1, x]

                if not xm and not ym:
                    largest_label = largest_label + 1
                    label_[y, x] = largest_label

                elif xm and not ym:
                    label_[y, x] = label_[y, x-1]

                elif not xm and ym:
                    label_[y, x] = label_[y-1, x]

                else:
                    union2(label_, y, label_[y-1, x], label_[y, x-1])
                    label_[y, x] = label_[y, x-1]


def hk2_2(lattice_, label_):
    largest_label = 0
    for y in range(ny):
        for x in range(nx):
            if lattice_[y, x] == 1:
                # Can I use other BD condition?
                if x == 0:
                    xm = 0
                else:
                    xm = lattice_[y, x-1]
                if y == 0:
                    ym = 0
                else:
                    ym = lattice_[y-1, x]

                last_connect = ''

                if not xm and not ym:
                    largest_label = largest_label + 1
                    label_[y, x] = largest_label

                else:
                    if xm:
                        label_[y, x] = label_[y, x-1]
                        last_connect = 'xm'

                    if ym:
                        if last_connect == '':
                            label_[y, x] = label_[y-1, x]
                        elif last_connect == 'xm':
                            union2(label_, y, label_[y-1, x], label_[y, x-1])
                            label_[y, x] = label_[y-1, x]


def union3(label_, z_, label1_, label2_):
    for iz in range(z_ + 1):
        for ix in range(nx):
            for iy in range(ny):
                if label_[iz, iy, ix] == label1_:
                    label_[iz, iy, ix] = label2_


def hk3(lattice_, label_):
    largest_label = 0
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                if lattice_[z, y, x] == 1:
                    # Can I use other BD condition?
                    xm = 0
                    if x > 0:
                        xm = lattice_[z, y, x-1]
                    ym = 0
                    if y > 0:
                        ym = lattice_[z, y-1, x]
                    zm = 0
                    if z > 0:
                        zm = lattice_[z-1, y, x]

                    if not xm and not ym and not zm:
                        largest_label = largest_label + 1
                        label_[z, y, x] = largest_label

                    elif xm and not ym and not zm:  # xm
                        label_[z, y, x] = label_[z, y, x-1]

                    elif not xm and ym and not zm:  # ym
                        label_[z, y, x] = label_[z, y-1, x]

                    elif not xm and not ym and zm:  # zm
                        label_[z, y, x] = label_[z-1, y, x]

                    elif not xm and ym and zm:  # ym and zm
                        union3(label_, z, label_[z-1, y, x], label_[z, y-1, x])
                        label_[z, y, x] = label_[z-1, y, x]

                    elif xm and not ym and zm:  # xm and zm
                        union3(label_, z, label_[z-1, y, x], label_[z, y, x-1])
                        label_[z, y, x] = label_[z-1, y, x]

                    elif xm and ym and not zm:  # xm and ym
                        union3(label_, z, label_[z, y-1, x], label_[z, y, x-1])
                        label_[z, y, x] = label_[z, y-1, x]

                    else:  # xm and ym and zm
                        union3(label_, z, label_[z-1, y, x], label_[z, y, x-1])
                        union3(label_, z, label_[z, y-1, x], label_[z, y, x-1])
                        label_[z, y, x] = label_[z-1, y, x]


def hk3_2(lattice_, label_):
    largest_label = 0
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                if lattice_[z, y, x] == 1:
                    # Can I use other BD condition?
                    xm = 0
                    if x > 0:
                        xm = lattice_[z, y, x - 1]
                    ym = 0
                    if y > 0:
                        ym = lattice_[z, y - 1, x]
                    zm = 0
                    if z > 0:
                        zm = lattice_[z - 1, y, x]

                    last_connect = ''

                    if not xm and not ym and not zm:
                        largest_label = largest_label + 1
                        label_[z, y, x] = largest_label

                    else:

                        if xm:  # handle xm first
                            label_[z, y, x] = label_[z, y, x-1]
                            last_connect = 'xm'

                        if ym:
                            if last_connect == '':
                                label_[z, y, x] = label_[z, y-1, x]
                            elif last_connect == 'xm':
                                union3(label_, z, label_[z, y-1, x], label_[z, y, x-1])
                                label_[z, y, x] = label_[z, y-1, x]
                            last_connect = 'ym'

                        if zm:
                            if last_connect == '':
                                label_[z, y, x] = label_[z-1, y, x]
                            elif last_connect == 'xm':
                                union3(label_, z, label_[z-1, y, x], label_[z, y, x-1])
                                label_[z, y, x] = label_[z-1, y, x]
                            elif last_connect == 'ym':
                                union3(label_, z, label_[z-1, y, x], label_[z, y-1, x])
                                label_[z, y, x] = label_[z-1, y, x]


def union4(label_, t_, label1_, label2_):
    for it in range(t_ + 1):
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    if label_[it, iz, iy, ix] == label1_:
                        label_[it, iz, iy, ix] = label2_


def hk4(lattice_, label_):
    largest_label = 0
    for t in range(nt):
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    if lattice_[t, z, y, x] == 1:
                        # Can I use other BD condition?
                        xm = 0
                        if x > 0:
                            xm = lattice_[t, z, y, x-1]
                        ym = 0
                        if y > 0:
                            ym = lattice_[t, z, y-1, x]
                        zm = 0
                        if z > 0:
                            zm = lattice_[t, z-1, y, x]
                        tm = 0
                        if t > 0:
                            tm = lattice_[t-1, z, y, x]

                        if not xm and not ym and not zm and not tm:
                            largest_label = largest_label + 1
                            label_[t, z, y, x] = largest_label

                        elif xm and not ym and not zm and not tm:  # xm
                            label_[t, z, y, x] = label_[t, z, y, x-1]

                        elif not xm and ym and not zm and not tm:  # ym
                            label_[t, z, y, x] = label_[t, z, y-1, x]

                        elif not xm and not ym and zm and not tm:  # zm
                            label_[t, z, y, x] = label_[t, z-1, y, x]

                        elif not xm and not ym and not zm and tm:  # tm
                            label_[t, z, y, x] = label_[t-1, z, y, x]

                        # ----------------------------------------------------------

                        elif xm and ym and not zm and not tm:  # xm and ym
                            union4(label_, t, label_[t, z, y-1, x], label_[t, z, y, x-1])
                            label_[t, z, y, x] = label_[t, z, y, x-1]

                        elif xm and not ym and zm and not tm:  # xm and zm
                            union4(label_, t, label_[t, z-1, y, x], label_[t, z, y, x-1])
                            label_[t, z, y, x] = label_[t, z, y, x-1]

                        elif not xm and ym and zm and not tm:  # ym and zm
                            union4(label_, t, label_[t, z-1, y, x], label_[t, z, y-1, x])
                            label_[t, z, y, x] = label_[t, z, y-1, x]

                        elif xm and not ym and not zm and tm:  # xm and tm
                            union4(label_, t, label_[t-1, z, y, x], label_[t, z, y, x-1])
                            label_[t, z, y, x] = label_[t, z, y, x-1]

                        elif not xm and ym and not zm and tm:  # ym and tm
                            union4(label_, t, label_[t-1, z, y, x], label_[t, z, y-1, x])
                            label_[t, z, y, x] = label_[t, z, y-1, x]

                        elif not xm and not ym and zm and tm:  # zm and tm
                            union4(label_, t, label_[t-1, z, y, x], label_[t, z-1, y, x])
                            label_[t, z, y, x] = label_[t, z-1, y, x]

                        # ----------------------------------------------------------

                        elif xm and ym and zm and not tm:  # xm and ym and zm
                            union4(label_, t, label_[t, z-1, y, x], label_[t, z, y, x-1])
                            union4(label_, t, label_[t, z, y-1, x], label_[t, z, y, x-1])
                            label_[t, z, y, x] = label_[t, z-1, y, x]

                        elif xm and ym and not zm and tm:  # xm and ym and tm
                            union4(label_, t, label_[t, z, y-1, x], label_[t, z, y, x-1])
                            union4(label_, t, label_[t-1, z, y, x], label_[t, z, y, x-1])
                            label_[t, z, y, x] = label_[t-1, z, y, x]

                        elif not xm and ym and zm and tm:  # ym and zm and tm
                            union4(label_, t, label_[t, z-1, y, x], label_[t, z, y-1, x])
                            union4(label_, t, label_[t-1, z, y, x], label_[t, z, y-1, x])
                            label_[t, z, y, x] = label_[t-1, z, y, x]

                        elif xm and not ym and zm and tm:  # xm and zm and tm
                            union4(label_, t, label_[t, z-1, y, x], label_[t, z, y, x-1])
                            union4(label_, t, label_[t-1, z, y, x], label_[t, z, y, x-1])
                            label_[t, z, y, x] = label_[t-1, z, y, x]
                        # ----------------------------------------------------------

                        else:  # xm and ym and zm and tm
                            union4(label_, t, label_[t, z, y-1, x], label_[t, z, y, x-1])
                            union4(label_, t, label_[t, z-1, y, x], label_[t, z, y, x-1])
                            union4(label_, t, label_[t-1, z, y, x], label_[t, z, y, x-1])
                            label_[t, z, y, x] = label_[t-1, z, y, x]


def hk4_2(lattice_, label_):
    largest_label = 0
    for t in range(nt):
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    if lattice_[t, z, y, x] == 1:
                        # Can I use other BD condition?
                        xm = 0
                        if x > 0:
                            xm = lattice_[t, z, y, x - 1]
                        ym = 0
                        if y > 0:
                            ym = lattice_[t, z, y - 1, x]
                        zm = 0
                        if z > 0:
                            zm = lattice_[t, z - 1, y, x]
                        tm = 0
                        if t > 0:
                            tm = lattice_[t - 1, z, y, x]

                        last_connect = ''

                        if not xm and not ym and not zm and not tm:
                            largest_label = largest_label + 1
                            label_[t, z, y, x] = largest_label

                        else:

                            if xm:  # handle xm first
                                label_[t, z, y, x] = label_[t, z, y, x-1]
                                last_connect = 'xm'

                            if ym:
                                if last_connect == '':
                                    label_[t, z, y, x] = label_[t, z, y-1, x]
                                elif last_connect == 'xm':
                                    union4(label_, t, label_[t, z, y-1, x], label_[t, z, y, x-1])
                                    label_[t, z, y, x] = label_[t, z, y-1, x]
                                last_connect = 'ym'

                            if zm:
                                if last_connect == '':
                                    label_[t, z, y, x] = label_[t, z-1, y, x]
                                elif last_connect == 'xm':
                                    union4(label_, t, label_[t, z-1, y, x], label_[t, z, y, x-1])
                                    label_[t, z, y, x] = label_[t, z-1, y, x]
                                elif last_connect == 'ym':
                                    union4(label_, t, label_[t, z-1, y, x], label_[t, z, y-1, x])
                                    label_[t, z, y, x] = label_[t, z-1, y, x]
                                last_connect = 'zm'

                            if tm:
                                if last_connect == '':
                                    label_[t, z, y, x] = label_[t-1, z, y, x]
                                elif last_connect == 'xm':
                                    union4(label_, t, label_[t-1, z, y, x], label_[t, z, y, x-1])
                                    label_[t, z, y, x] = label_[t-1, z, y, x]
                                elif last_connect == 'ym':
                                    union4(label_, t, label_[t-1, z, y, x], label_[t, z, y-1, x])
                                    label_[t, z, y, x] = label_[t-1, z, y, x]
                                elif last_connect == 'zm':
                                    union4(label_, t, label_[t-1, z, y, x], label_[t, z-1, y, x])
                                    label_[t, z, y, x] = label_[t-1, z, y, x]


lattice = np.random.randint(0, 2, nx*ny).reshape(ny, nx)
label = np.random.randint(0, 1, nx*ny).reshape(ny, nx)

hk2(lattice, label)
print(lattice)
print(label)
print('...')

bak = label.copy()
hk2_2(lattice, label)
print(label)
print('...')

print(bak - label)
print(np.sum(bak - label))
print('...')

lattice = np.random.randint(0, 2, nx*ny*nz).reshape((nz, ny, nx))
lattice[..., 1:] = 0
label = np.random.randint(0, 1, nx*ny*nz).reshape((nz, ny, nx))

hk3(lattice, label)
print(lattice[..., 0])
print(label[..., 0])
print('...')
bak = label.copy()

hk3_2(lattice, label)
print(lattice[..., 0])
print(label[..., 0])
print('...')

print((label-bak)[..., 0])
print('...')

lattice = np.random.randint(0, 2, nx*ny*nz).reshape((nz, ny, nx))
hk3(lattice, label)
bak = label.copy()
hk3_2(lattice, label)
print(np.sum(bak - label))
print('...')

random = np.random.randint(0, 2, nx*ny*nz*nt).reshape((nt, nz, ny, nx))
lattice = np.zeros(shape=(nt, nz, ny, nx),dtype='int')
lattice[..., 0, 0] = random[..., 0, 0]
label = lattice * 0

hk4(lattice, label)
print(lattice[..., 0, 0])
print(label[..., 0, 0])
print('...')
bak = label.copy()

hk4_2(lattice, label)
print(lattice[..., 0, 0])
print(label[..., 0, 0])
print('...')

print((label-bak)[..., 0, 0])

lattice = np.random.randint(0, 2, nx*ny*nz*nt).reshape((nt, nz, ny, nx))

hk4(lattice, label)
bak = label.copy()
hk4_2(lattice, label)
print(np.sum(bak - label))
