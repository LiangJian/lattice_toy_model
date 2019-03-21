import numba
from dim4_pos import *


@numba.jit
def union4(label_, t_, label1_, label2_):
    for it_ in range(t_ + 1):
        for iz_ in range(nz):
            for iy_ in range(ny):
                for ix_ in range(nx):
                    if label_[it_, iz_, iy_, ix_] == label1_:
                        label_[it_, iz_, iy_, ix_] = label2_


@numba.jit
def hk4_prop(lattice_, label_, value_, props_):
    largest_label = 0
    for t in range(nt):
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    if lattice_[t, z, y, x] == value_:
                        # Can I use other BD condition?
                        xm = False
                        if x > 0:
                            xm = lattice_[t, z, y, x - 1] == value_ and props_[t, z, y, x-1, 3]
                        ym = False
                        if y > 0:
                            ym = lattice_[t, z, y - 1, x] == value_ and props_[t, z, y-1, x, 2]
                        zm = False
                        if z > 0:
                            zm = lattice_[t, z - 1, y, x] == value_ and props_[t, z-1, y, x, 1]
                        tm = False
                        if t > 0:
                            tm = lattice_[t - 1, z, y, x] == value_ and props_[t-1, z, y, x, 0]

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


@numba.jit
def hk4_2_prop(lattice_, label_, value_, props_):
    largest_label = 0
    for t in range(nt):
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    if lattice_[t, z, y, x] == value_:
                        # Can I use other BD condition?
                        xm = False
                        if x > 0:
                            xm = lattice_[t, z, y, x - 1] == value_ and props_[t, z, y, x-1, 3]
                        ym = False
                        if y > 0:
                            ym = lattice_[t, z, y - 1, x] == value_ and props_[t, z, y-1, x, 2]
                        zm = False
                        if z > 0:
                            zm = lattice_[t, z - 1, y, x] == value_ and props_[t, z-1, y, x, 1]
                        tm = False
                        if t > 0:
                            tm = lattice_[t - 1, z, y, x] == value_ and props_[t-1, z, y, x, 0]

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


@numba.jit
def handle_labels(labels_):
    """ This function is only for test use
    """
    label_list_ = []
    cluster_points_ = {}
    for it_ in range(nt):
        for iz_ in range(nz):
            for iy_ in range(ny):
                for ix_ in range(nx):
                    label_ = int(labels_[it_, iz_, iy_, ix_])
                    if label_ > 0:
                        if label_ not in label_list_:
                            label_list_.append(label_)
                            cluster_points_[label_] = []
                            cluster_points_[label_].append((it_, iz_, iy_, ix_))
                        else:
                            cluster_points_[label_].append((it_, iz_, iy_, ix_))
    return np.array(label_list_), cluster_points_
