import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def data_coord2view_coord(data, reso, data_min, data_max):
    data_range = data_max - data_min
    plt.plot((data - data_min) / data_range * reso)
    plt.show()
    
    return (data - data_min) / data_range * reso


def nearest_neighbours_1d(xs, ys, reso, n_neighbours):
    im = np.zeros([reso])
    extent = [np.min(xs), np.max(xs)]
    
    print(im)
    print(extent)
    print()
    print()

    xv = data_coord2view_coord(xs, reso, extent[0], extent[1])
    for x in range(reso):
        xp = (xv - x)

        d = xp

        im[x] = 1 / np.sum(d[np.argpartition(d.ravel(), n_neighbours)[:n_neighbours]])
    im = im[np.newaxis, :]
    print(im)
    return im, extent

def nearest_neighbours(xs, ys, reso, n_neighbours):
    im = np.zeros([reso, reso])
    extent = [np.min(xs), np.max(xs), np.min(ys), np.max(ys)]

    xv = data_coord2view_coord(xs, reso, extent[0], extent[1])
    yv = data_coord2view_coord(ys, reso, extent[2], extent[3])
    for x in range(reso):
        for y in range(reso):
            xp = (xv - x)
            yp = (yv - y)

            d = np.sqrt(xp**2 + yp**2)
            print(d.shape)
            print(d.ravel().shape)

            exit()
            print(d[np.argpartition(d.ravel(), n_neighbours)[:n_neighbours]].shape)
            im[y][x] = 1 / np.sum(d[np.argpartition(d.ravel(), n_neighbours)[:n_neighbours]])

    return im, extent

n = 1000
xs = np.random.randn(n)
ys = np.random.randn(n)
resolution = 250

fig, axes = plt.subplots(2, 2)

for ax, neighbours in zip(axes.flatten(), [0, 16]):
    if neighbours == 0:
        ax.plot(xs, ys, 'k.', markersize=2)
        ax.set_aspect('equal')
        ax.set_title("Scatter Plot")
    else:
        im, extent = nearest_neighbours_1d(xs, ys, resolution, neighbours)
        print(im)
        ax.imshow(im,  cmap=cm.jet)
        ax.set_title("Smoothing over %d neighbours" % neighbours)
        # ax.set_xlim(extent[0], extent[1])
        # ax.set_ylim(0,1)
plt.show()