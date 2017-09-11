import numpy as np

'''
To test Hoshen–Kopelman algorithm
'''

nx = 8
ny = 8
lattice = np.random.randint(0, 2, nx*ny).reshape(nx, ny)
label = np.random.randint(0, 1, nx*ny).reshape(nx, ny)


def union(y_, label1_, label2_):
    for ix in range(nx):
        for iy in range(y_ + 1):
            if label[ix, iy] == label1_:
                label[ix, iy] = label2_


largest_label = 0
for y in range(ny):
    for x in range(nx):
        if lattice[x, y]:
            if x == 0:
                left = 0
            else:
                left = lattice[x-1, y]
            if y == 0:
                above = 0
            else:
                above = lattice[x, y-1]
            if not left and not above:
                largest_label = largest_label + 1
                label[x, y] = largest_label
            elif left and not above:
                label[x, y] = label[x-1, y]
            elif not left and above:
                label[x, y] = label[x, y-1]
            else:
                union(y, label[x, y-1], label[x-1, y])
                label[x, y] = label[x-1, y]

print(lattice)
print(label)
