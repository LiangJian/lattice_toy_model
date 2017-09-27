#######################################
# Hoshen–Kopelman algorithm for 4 dim
#######################################


from dim4_pos import *
import numpy as np


def equive(pos_, labels_, label1_, label2_):
    for it_ in range(pos_[0] + 1):
        for iz_ in range(nz):
            for iy_ in range(ny):
                for ix_ in range(nx):
                    if labels_[it_, iz_, iy_, ix_] == label1_:
                        labels_[it_, iz_, iy_, ix_] = label2_


def bond(field_, value_):
    labels_ = np.zeros(shape=field_.shape)
    largest_label = 0
    for ix_ in sites:
        index_ = index[ix_]
        if field_[tuple(index_)] == value_:  # choose 1 or -1
            xm = np.zeros(4, dtype='int')
            for iu in range(4):
                if index_[iu] == 0:  # bound
                    xm[iu] = 0
                else:
                    index_[iu] -= 1
                    if field_[tuple(index_)] == value_:
                        xm[iu] = 1
                    index_[iu] += 1

                # no bond
                if not xm[0] and not xm[1] and not xm[2] and not xm[3]:
                    largest_label = largest_label + 1
                    labels_[tuple(index_)] = largest_label

                # one bond
                # t
                elif xm[0] and not xm[1] and not xm[2] and not xm[3]:
                    labels_[tuple(index_)] = labels_[index_[0] - 1, index_[1], index_[2], index_[3]]
                # z
                elif not xm[0] and xm[1] and not xm[2] and not xm[3]:
                    labels_[tuple(index_)] = labels_[index_[0], index_[1] - 1, index_[2], index_[3]]
                # y
                elif not xm[0] and not xm[1] and xm[2] and not xm[3]:
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2] - 1, index_[3]]
                # x
                elif not xm[0] and not xm[1] and not xm[2] and xm[3]:
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2], index_[3] - 1]

                # two bonds
                # Since only the outermost loop may be scanned by part, order t z y x is somehow better
                # t z
                elif xm[0] and xm[1] and not xm[2] and not xm[3]:
                    equive(tuple(index_), labels_,
                           labels_[index_[0] - 1, index_[1], index_[2], index_[3]],
                           labels_[index_[0], index_[1] - 1, index_[2], index_[3]])
                    labels_[tuple(index_)] = labels_[index_[0], index_[1] - 1, index_[2], index_[3]]

                # t y
                elif xm[0] and not xm[1] and xm[2] and not xm[3]:
                    equive(tuple(index_), labels_,
                           labels_[index_[0] - 1, index_[1], index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2] - 1, index_[3]])
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2] - 1, index_[3]]
                # t x
                elif xm[0] and not xm[1] and not xm[2] and xm[3]:
                    equive(tuple(index_), labels_,
                           labels_[index_[0] - 1, index_[1], index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3]] - 1)
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2], index_[3]] - 1
                # z y
                elif not xm[0] and xm[1] and xm[2] and not xm[3]:
                    equive(tuple(index_), labels_,
                           labels_[index_[0], index_[1] - 1, index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2] - 1, index_[3]])
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2] - 1, index_[3]]
                # z x
                elif not xm[0] and xm[1] and not xm[2] and xm[3]:
                    equive(tuple(index_), labels_,
                           labels_[index_[0], index_[1] - 1, index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2], index_[3] - 1]
                # y x
                elif not xm[0] and not xm[1] and xm[2] and xm[3]:
                    equive(tuple(index_), labels_,
                           labels_[index_[0], index_[1], index_[2] - 1, index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2], index_[3] - 1]

                # three bonds
                # t z y
                elif xm[0] and xm[1] and xm[2] and not xm[3]:
                    equive(tuple(index_), labels_,
                           labels_[index_[0] - 1, index_[1], index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2] - 1, index_[3]])
                    equive(tuple(index_), labels_,
                           labels_[index_[0], index_[1] - 1, index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2] - 1, index_[3]])
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2] - 1, index_[3]]
                # t z x
                elif xm[0] and xm[1] and not xm[2] and xm[3]:
                    equive(tuple(index_), labels_,
                           labels_[index_[0] - 1, index_[1], index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    equive(tuple(index_), labels_,
                           labels_[index_[0], index_[1] - 1, index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2], index_[3] - 1]
                # t y x
                elif xm[0] and not xm[1] and xm[2] and xm[3]:
                    equive(tuple(index_), labels_,
                           labels_[index_[0] - 1, index_[1], index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    equive(tuple(index_), labels_,
                           labels_[index_[0], index_[1], index_[2] - 1, index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2], index_[3] - 1]
                # z y x
                elif not xm[0] and xm[1] and xm[2] and xm[3]:
                    equive(tuple(index_), labels_,
                           labels_[index_[0], index_[1] - 1, index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    equive(tuple(index_), labels_,
                           labels_[index_[0], index_[1], index_[2] - 1, index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2], index_[3] - 1]

                # four bonds
                # t z y x
                else:
                    equive(tuple(index_), labels_,
                           labels_[index_[0] - 1, index_[1], index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    equive(tuple(index_), labels_,
                           labels_[index_[0], index_[1] - 1, index_[2], index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    equive(tuple(index_), labels_,
                           labels_[index_[0], index_[1], index_[2] - 1, index_[3]],
                           labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                    labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2], index_[3] - 1]
    return labels_


def handle_labels(labels_):
    """ This function is only for test use
    """
    label_list = []
    cluster_points = {}
    for it_ in range(nt):
        for iz_ in range(nz):
            for iy_ in range(ny):
                for ix_ in range(nx):
                    label_ = int(labels_[it_, iz_, iy_, ix_])
                    if label_ > 0:
                        if label_ not in label_list:
                            label_list.append(label_)
                            cluster_points[label_] = []
                        else:
                            cluster_points[label_].append((it_, iz_, iy_, ix_))
    return label_list, cluster_points


def count(labels_):
    """ This function is only for test use
    """
    count_ = 0
    for it_ in range(nt):
        for iz_ in range(nz):
            for iy_ in range(ny):
                for ix_ in range(nx):
                    # print(labels_[it_, iz_, iy_, ix_])
                    if labels_[it_, iz_, iy_, ix_] > 0:
                        count_ += 1
    return count_
