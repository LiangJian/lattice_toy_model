#######################################
# Hoshen–Kopelman algorithm for 4 dim
#######################################


from dim4_pos import *
import numpy as np
import copy


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
        index_ = list(index[ix_])
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
                       labels_[index_[0], index_[1], index_[2], index_[3] - 1])
                labels_[tuple(index_)] = labels_[index_[0], index_[1], index_[2], index_[3] - 1]
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


def handle_one_bond(x_, props_, largest_label, labels_, index_):
    index0_ = copy.deepcopy(index_)
    index0_[x_] -= 1
    if np.random.random() < props_[tuple(index0_)][x_]:
        labels_[tuple(index_)] = labels_[tuple(index0_)]
    else:
        largest_label[0] += 1
        labels_[tuple(index_)] = largest_label[0]


def handle_two_bonds(x_, y_, props_, largest_label, labels_, index_):
    p0_ = np.random.random()
    p1_ = np.random.random()
    index0_ = copy.deepcopy(index_)
    index1_ = copy.deepcopy(index_)
    index0_[x_] -= 1
    index1_[y_] -= 1
    if p0_ < props_[tuple(index0_)][x_] and p1_ < props_[tuple(index1_)][y_]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index1_)])
        labels_[tuple(index_)] = labels_[tuple(index1_)]
    elif p0_ < props_[tuple(index0_)][x_]:
        labels_[tuple(index_)] = labels_[tuple(index0_)]
    elif p1_ < props_[tuple(index1_)][y_]:
        labels_[tuple(index_)] = labels_[tuple(index1_)]
    else:
        largest_label[0] += 1
        labels_[tuple(index_)] = largest_label[0]


def handle_three_bonds(x_, y_, z_, props_, largest_label, labels_, index_):
    p0_ = np.random.random()
    p1_ = np.random.random()
    p2_ = np.random.random()
    index0_ = copy.deepcopy(index_)
    index1_ = copy.deepcopy(index_)
    index2_ = copy.deepcopy(index_)
    index0_[x_] -= 1
    index1_[y_] -= 1
    index2_[z_] -= 1
    if p0_ < props_[tuple(index0_)][x_] and p1_ < props_[tuple(index1_)][y_] and p2_ < props_[tuple(index2_)][z_]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index1_)])
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index2_)])
        labels_[tuple(index_)] = labels_[tuple(index2_)]
    elif p0_ < props_[tuple(index0_)][x_] and p1_ < props_[tuple(index1_)][y_]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index1_)])
        labels_[tuple(index_)] = labels_[tuple(index1_)]
    elif p0_ < props_[tuple(index0_)][x_] and p2_ < props_[tuple(index2_)][z_]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index2_)])
        labels_[tuple(index_)] = labels_[tuple(index2_)]
    elif p1_ < props_[tuple(index1_)][y_] and p2_ < props_[tuple(index2_)][z_]:
        equive(tuple(index_), labels_, labels_[tuple(index1_)], labels_[tuple(index2_)])
        labels_[tuple(index_)] = labels_[tuple(index2_)]
    elif p0_ < props_[tuple(index0_)][x_]:
        labels_[tuple(index_)] = labels_[tuple(index0_)]
    elif p1_ < props_[tuple(index1_)][y_]:
        labels_[tuple(index_)] = labels_[tuple(index1_)]
    elif p2_ < props_[tuple(index2_)][z_]:
        labels_[tuple(index_)] = labels_[tuple(index2_)]
    else:
        largest_label[0] += 1
        labels_[tuple(index_)] = largest_label[0]


def handle_four_bonds(props_, largest_label, labels_, index_):
    p0_ = np.random.random()
    p1_ = np.random.random()
    p2_ = np.random.random()
    p3_ = np.random.random()
    index0_ = copy.deepcopy(index_)
    index1_ = copy.deepcopy(index_)
    index2_ = copy.deepcopy(index_)
    index3_ = copy.deepcopy(index_)
    index0_[0] -= 1
    index1_[1] -= 1
    index2_[2] -= 1
    index3_[3] -= 1
    if p0_ < props_[tuple(index0_)][0] and p1_ < props_[tuple(index1_)][1] and p2_ < props_[tuple(index2_)][2] and \
            p3_ < props_[tuple(index3_)][3]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index1_)])
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index2_)])
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index3_)])
        labels_[tuple(index_)] = labels_[tuple(index3_)]
    elif p0_ < props_[tuple(index0_)][0] and p1_ < props_[tuple(index1_)][1] and p2_ < props_[tuple(index2_)][2]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index1_)])
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index2_)])
        labels_[tuple(index_)] = labels_[tuple(index2_)]
    elif p0_ < props_[tuple(index0_)][0] and p1_ < props_[tuple(index1_)][1] and p3_ < props_[tuple(index3_)][3]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index1_)])
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index3_)])
        labels_[tuple(index_)] = labels_[tuple(index3_)]
    elif p0_ < props_[tuple(index0_)][0] and p2_ < props_[tuple(index2_)][2] and p3_ < props_[tuple(index3_)][3]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index2_)])
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index3_)])
        labels_[tuple(index_)] = labels_[tuple(index3_)]
    elif p1_ < props_[tuple(index1_)][1] and p2_ < props_[tuple(index2_)][2] and p3_ < props_[tuple(index3_)][3]:
        equive(tuple(index_), labels_, labels_[tuple(index1_)], labels_[tuple(index2_)])
        equive(tuple(index_), labels_, labels_[tuple(index1_)], labels_[tuple(index3_)])
        labels_[tuple(index_)] = labels_[tuple(index3_)]
    elif p0_ < props_[tuple(index0_)][0] and p1_ < props_[tuple(index1_)][1]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index1_)])
        labels_[tuple(index_)] = labels_[tuple(index1_)]
    elif p0_ < props_[tuple(index0_)][0] and p2_ < props_[tuple(index2_)][2]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index2_)])
        labels_[tuple(index_)] = labels_[tuple(index2_)]
    elif p0_ < props_[tuple(index0_)][0] and p3_ < props_[tuple(index3_)][3]:
        equive(tuple(index_), labels_, labels_[tuple(index0_)], labels_[tuple(index3_)])
        labels_[tuple(index_)] = labels_[tuple(index3_)]
    elif p1_ < props_[tuple(index1_)][1] and p2_ < props_[tuple(index2_)][2]:
        equive(tuple(index_), labels_, labels_[tuple(index1_)], labels_[tuple(index2_)])
        labels_[tuple(index_)] = labels_[tuple(index2_)]
    elif p1_ < props_[tuple(index1_)][1] and p3_ < props_[tuple(index3_)][3]:
        equive(tuple(index_), labels_, labels_[tuple(index1_)], labels_[tuple(index3_)])
        labels_[tuple(index_)] = labels_[tuple(index3_)]
    elif p2_ < props_[tuple(index2_)][2] and p3_ < props_[tuple(index3_)][3]:
        equive(tuple(index_), labels_, labels_[tuple(index2_)], labels_[tuple(index3_)])
        labels_[tuple(index_)] = labels_[tuple(index3_)]
    elif p0_ < props_[tuple(index0_)][0]:
        labels_[tuple(index_)] = labels_[tuple(index0_)]
    elif p1_ < props_[tuple(index1_)][1]:
        labels_[tuple(index_)] = labels_[tuple(index1_)]
    elif p2_ < props_[tuple(index2_)][2]:
        labels_[tuple(index_)] = labels_[tuple(index2_)]
    elif p3_ < props_[tuple(index3_)][3]:
        labels_[tuple(index_)] = labels_[tuple(index3_)]
    else:
        largest_label[0] += 1
        labels_[tuple(index_)] = largest_label[0]


def bond_prop(field_, value_, props_):
    labels_ = np.zeros(shape=field_.shape)
    largest_label = 0
    for ix_ in sites:
        index_ = list(index[ix_])
        if field_[tuple(index_)] == value_:  # choose 1 or -1
            xm = np.zeros(4, dtype='int')
            for iu in range(4):
                if index_[iu] == 0:
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
                largest_label_tmp_list = [largest_label]
                handle_one_bond(0, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # z
            elif not xm[0] and xm[1] and not xm[2] and not xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_one_bond(1, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # y
            elif not xm[0] and not xm[1] and xm[2] and not xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_one_bond(2, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # x
            elif not xm[0] and not xm[1] and not xm[2] and xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_one_bond(3, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # two bonds
            # Since only the outermost loop may be scanned by part, order t z y x is somehow better
            # t z
            elif xm[0] and xm[1] and not xm[2] and not xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_two_bonds(0, 1, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # t y
            elif xm[0] and not xm[1] and xm[2] and not xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_two_bonds(0, 2, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # t x
            elif xm[0] and not xm[1] and not xm[2] and xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_two_bonds(0, 3, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # z y
            elif not xm[0] and xm[1] and xm[2] and not xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_two_bonds(1, 2, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # z x
            elif not xm[0] and xm[1] and not xm[2] and xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_two_bonds(1, 3, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # y x
            elif not xm[0] and not xm[1] and xm[2] and xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_two_bonds(2, 3, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # three bonds
            # t z y
            elif xm[0] and xm[1] and xm[2] and not xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_three_bonds(0, 1, 2, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # t z x
            elif xm[0] and xm[1] and not xm[2] and xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_three_bonds(0, 1, 3, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # t y x
            elif xm[0] and not xm[1] and xm[2] and xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_three_bonds(0, 2, 3, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # z y x
            elif not xm[0] and xm[1] and xm[2] and xm[3]:
                largest_label_tmp_list = [largest_label]
                handle_three_bonds(1, 2, 3, props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

            # four bonds
            # t z y x
            else:
                largest_label_tmp_list = [largest_label]
                handle_four_bonds(props_, largest_label_tmp_list, labels_, index_)
                largest_label = largest_label_tmp_list[0]

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
                            cluster_points[label_].append((it_, iz_, iy_, ix_))
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
                    if labels_[it_, iz_, iy_, ix_] > 0:
                        count_ += 1
    return count_
