#######################################
# Hoshen–Kopelman algorithm for 4 dim
# new version
#######################################


from dim4_pos import *
import numpy as np


def handle_connections(index_, xm, labels_list_, labels_, largest_label):
    if not (xm[0] or xm[1] or xm[2] or xm[3]):
        labels_list_[largest_label[0]] = []
        labels_list_[largest_label[0]].append(index_)
        labels_[index_] = largest_label[0]
        largest_label[0] = largest_label[0] + 1
    else:
        if xm[0]:
            labels_list_[labels_[move_back(index_, 0)]] += [index_]
            labels_[index_] = labels_[move_back(index_, 0)]
            if xm[1]:
                pass




def bond_prop_new(field_, value_, props_):
    labels_list_ = {}
    labels_ = np.zeros(shape=field_.shape, dtype='int')
    largest_label = [0]
    for ix_ in sites:
        index_ = index[ix_]
        if field_[index_] == value_:  # choose 1 or -1. Do I really need this?
            xm = np.array([False] * 4)  # if there is connection at x_minus
            for iu_ in range(4):
                if index_[iu_] != 0:  # not reach the boundary
                    if field_[move_back(index_, iu_)] == value_:  # backward elements carries the same value
                        if np.random.random() < props_[move_back(index_, iu_)][iu_]:   # probability of the connection
                            xm[iu_] = True  #
            handle_connections(index_, xm, labels_list_, labels_, largest_label)
    return labels_, labels_list_
