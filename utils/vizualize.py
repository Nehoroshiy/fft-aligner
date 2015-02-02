from matplotlib import rc
import numpy as np
import matplotlib.pyplot as plt

#Todo refac this - bad codex
def viz_way2D(seq1, seq2, alignment_seq1, alignment_seq2, way_mtx, way):
    way_show_mtx = np.zeros(shape=way_mtx.shape)
    for i in range(way_show_mtx.shape[0]):
        for j in range(way_show_mtx.shape[1]):
            way_show_mtx[i][j] = 0.7
    for i, j in way: way_show_mtx[i, j] = 0.8

    way_show_mtx[0, way_show_mtx.shape[1]-1] = 0
    way_show_mtx[way_show_mtx.shape[0]-1, 0] = 1

    plt.figure()
    plt.imshow(np.array(way_show_mtx), interpolation='nearest')

    width = len(way_show_mtx)
    height = len(way_show_mtx[0])

    for x in xrange(width):
        for y in xrange(height):
            plt.annotate(str(way_mtx[x][y]), xy=(y, x),  horizontalalignment='center', verticalalignment='center')

    plt.yticks(range(len(seq1)+1), ' ' + seq1)
    plt.xticks(range(len(seq2)+1), ' ' + seq2)

    font = {'family': 'monospace', 'size': 14}

    rc('font', **font)

    plt.title('%s\n%s' % (alignment_seq1, alignment_seq2))

    plt.show()