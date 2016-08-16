import numpy as np
import pickle
from configs import Locations
from itertools import product, combinations
from string import ascii_uppercase
# TODO: remove the reference to chiffatools
from chiffatools.linalg_routines import show_matrix_with_names, hierchical_clustering
from pre_processer import group_plot
from copy import deepcopy
from matplotlib import pyplot as plt
from scipy import ndimage

current_stack = pickle.load(open(Locations['dump']['raw']))
growth_speed, peak_growth = pickle.load(open(Locations['dump']['times']))
pad = pickle.load(open(Locations['pad']))

# TODO: refactor to disentangle namespaces


def compare_growth_curves(growth_curve1, growth_curve2):
    return np.sum(np.abs(growth_curve1 - growth_curve2))


def localize_strain():
    integer_unfold = np.zeros((60, 2)).astype(np.intp)
    flat_labels = np.zeros((60, )).astype(np.str)
    integer_fold = np.zeros((8, 12))
    flat_names = np.zeros((60, )).astype(np.str)

    for i, j in product(range(1, 7), range(1, 11)):
        flat_number = (i - 1) * 10 + j - 1
        integer_unfold[flat_number, :] = i, j
        integer_fold[i, j] = flat_number
        flat_labels[flat_number] = np.array('%s, %s : %s' % (ascii_uppercase[i], j + 1, pad[i, j]))
        flat_names[flat_number] = pad[i, j]

    return integer_fold, integer_unfold, flat_labels, flat_names


def plot_and_group(padded_table, group_pad=False):
    _, _, flat_labels, flat_names = localize_strain()
    show_matrix_with_names(padded_table, flat_labels.tolist(), flat_labels.tolist())
    if not group_pad:
        index = hierchical_clustering(padded_table, flat_names.tolist())
    else:
        index = np.argsort(flat_names)
    padded_table = padded_table[index, :]
    padded_table = padded_table[:, index]
    flat_labels = flat_labels[index]
    show_matrix_with_names(padded_table, flat_labels.tolist(), flat_labels.tolist())


def l2_dist(plates_stack, group_pad=False):
    padded_table = np.zeros((60, 60))

    _, integer_unfold, _, _ = localize_strain()

    for i, j in combinations(range(0, 60), 2):
        i1, j1 = integer_unfold[i]
        i2, j2 = integer_unfold[j]
        data1 = plates_stack[:, i1, j1]
        data2 = plates_stack[:, i2, j2]
        c_growth = compare_growth_curves(data1, data2)
        padded_table[i, j] = c_growth
        padded_table[j, i] = c_growth

    plot_and_group(padded_table, group_pad)


def l2_growth_difference(growth_speed_matrix, group_pad=False):
    padded_table = np.zeros((60, 60))

    _, integer_unfold, _, _ = localize_strain()

    for i, j in combinations(range(0, 60), 2):
        i1, j1 = integer_unfold[i]
        i2, j2 = integer_unfold[j]
        data1 = growth_speed_matrix[i1, j1]
        data2 = growth_speed_matrix[i2, j2]
        c_growth = np.abs(data1 - data2)
        padded_table[i, j] = c_growth
        padded_table[j, i] = c_growth

    plot_and_group(padded_table,  group_pad)


def group_by(group_ranges=(range(1, 10), range(10, 16))):

    integer_fold, integer_unfold, flat_labels, flat_names = localize_strain()

    repetability_score = {}
    unique_names = np.unique(flat_names)
    for name in unique_names:
        repetability_score[name] = [0, 0]
        acc = 0
        for i, j in combinations(np.nonzero(flat_names == name)[0], 2):
            i1, j1 = integer_unfold[i]
            i2, j2 = integer_unfold[j]
            data1 = current_stack[:, i1, j1]
            data2 = current_stack[:, i2, j2]
            repetability_score[name][0] += compare_growth_curves(data1, data2)
            data1 = growth_speed[i1, j1]
            data2 = growth_speed[i2, j2]
            repetability_score[name][1] += np.abs(data1 - data2)/np.mean(np.array(data1, data2))*100
            acc += 1
        repetability_score[name][0] /= acc
        repetability_score[name][1] /= acc

    for key, (val1, val2) in sorted(repetability_score.iteritems(), key=lambda x: x[1][1]):
        print "group %s \t  %s shape \t %s%% division speed" % (key, '{0:.2f}'.format(val1), '{0:.2f}'.format(val2))

    int_dict = {}
    for name in unique_names:
        is_int = True
        try:
            float(name)
        except ValueError:
            is_int = False
        if is_int:
            int_dict[int(float(name))] = name

    collector = np.zeros((len(group_ranges),)).tolist()
    collector = [[] for _ in collector]
    zoomlist = deepcopy(collector)
    for key, val in int_dict.iteritems():
        for i, _range in enumerate(group_ranges):
            if key in _range:
                zoomlist[i] = zoomlist[i] + np.nonzero(flat_names == val)[0].tolist()
                collector[i].append(repetability_score[val])
    collector = [np.mean(np.array(collection), axis=0).tolist() for collection in collector]
    for i, item in enumerate(collector):
        print i, "%s shape \t %s%% division speed" % ('{0:.2f}'.format(item[0]), '{0:.2f}'.format(item[1]))
    item = repetability_score['WT-ref']
    print 'WT-ref', "%s shape \t %s%% division speed" % ('{0:.2f}'.format(item[0]), '{0:.2f}'.format(item[1]))

    zoomlist.append(np.nonzero(flat_names == 'WT-ref')[0])
    zoomlist.append([1, 2])  # TODO: this is an illustration, we need to remove it in the future.
    zoomlist = [[tuple(integer_unfold[elt]) for elt in lst] for lst in zoomlist]

    group_plot(current_stack, zoomlist)


def regress_variability(division_speed, name_pad):

    integer_fold, _, _, _ = localize_strain()

    footprint = np.array([[0, 1, 0],
                          [1, 1, 1],
                          [0, 1, 0]])

    filter_index = name_pad != 'WT-ref'
    division_speed[filter_index] = np.NaN
    division_speed = ndimage.generic_filter(division_speed, np.nanmean, footprint=footprint)
    lvl2_filter = integer_fold == 0
    lvl2_filter[1, 1] = False
    division_speed[lvl2_filter] = np.nan
    division_speed = (division_speed / np.median(division_speed) - 1)*100
    plt.imshow(division_speed, interpolation='nearest')
    plt.colorbar()
    plt.savefig(Locations['output']['plate'], dpi=300)
    plt.show()


if __name__ == "__main__":
    # l2_dist(current_stack)
    # l2_growth_difference(growth_speed)
    # l2_dist(current_stack, group_pad=True)
    # l2_growth_difference(growth_speed, group_pad=False)
    # l2_growth_difference(growth_speed, group_pad=True)
    group_by()
    # regress_variability(growth_speed, pad)