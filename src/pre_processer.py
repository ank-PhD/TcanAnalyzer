__author__ = 'ank'

import os
import xlrd
import numpy as np
from itertools import product
from matplotlib import pyplot as plt
import matplotlib as mlb
from scipy.ndimage.filters import gaussian_filter1d as gf_1d
from collections import defaultdict
from time import time
from string import ascii_lowercase
from sklearn.gaussian_process import GaussianProcess
import pickle
from configs import Locations, base_folder_name
from src.CalibrateTCAN import load_corrector

# TODO: use cummulative growth speed to prevent too quick of the regression
# TODO: incorporate pad plotting

# TODO: Calculate the lag time as the difference in time that it took to reach the fastest growing point and the time it
#    would have taken if the yeast was growing from the initial OD to the OD at the fastest point at the top speed.

# TODO: import a exponential fitting module from the analysis that was used in the Jin's historical data

# TODO: debug the bump appearance in the early stages => DONE; the problem is in the original data.
# SLOPE: use more than 6 points used to compute the slope => Problem; it over-smoothes it and changes the precision of slope estimate

# TODO: create an OD block, where the OD variation is too likely to lead to a perturbed division speed, exclude it from the possible location where the OD gradient can be calculated

debug = False

tinit = time()
mlb.rcParams['font.size'] = 10.0
mlb.rcParams['figure.figsize'] = (30, 20)

file_location = 'U:/ank/2015/TcanScreen/04.27.2015'
file_name = 'Book1.xlsx'
pad_location = 'U:/ank/2015/TcanScreen/03.26.2015'
pad_name = 'pad.xlsx'
d_time = 15./60.

time_map = defaultdict(int)

def extract_plate(a1_coordinates, sheet, string=False):
    plate = np.zeros((8, 12))
    if string:
        plate = plate.astype(np.str)
    for i, j in product(range(0, 8), range(0, 12)):
        _i, _j = (np.array([i, j]) + np.array(list(a1_coordinates))).astype(np.uint32).tolist()
        # print _i, _j, sheet.cell(_i, _j).value
        if sheet.cell(_i, _j).value:
            plate[i, j] = sheet.cell(_i, _j).value
        else:
            print 'missing data @', _i, _j
            plate[i, j] = -1
    return plate


def extract_plate_dit(current_path):
    wb = xlrd.open_workbook(current_path)
    collect_dict = []
    for s in wb.sheets():
        if 'Magellan' in s.name:
            latest = 0
            for row in range(1, s.nrows):
                if row % 9 == 1:
                    plate = extract_plate((row, 1), s)
                    collect_dict.append(plate)
                    latest += d_time
    return np.array(collect_dict)


def read_pad(pad_path):
    print pad_path
    wb = xlrd.open_workbook(pad_path)
    plate = []
    for s in wb.sheets():
        if 'Sheet1' in s.name:
            plate = extract_plate((0, 0), s, string=True)
    pickle.dump(plate, open(Locations['pad'], 'w'))
    return plate


def smooth_plate(plate, window):
    re_plate = plate.copy()
    for i, j in product(range(0, 8), range(0, 12)):
        re_plate[:, i, j] = gf_1d(plate[:, i, j], window)
    return re_plate


def plot_growth(plates_stack, grad=False, dumpType=None, NoShow=False):
    fig = plt.gcf()
    fig.canvas.set_window_title('min:%s, max:%s, total points: %s' % (np.min(plates_stack),
                                                                      np.max(plates_stack),
                                                                      plates_stack.shape[0]))
    speeds = np.zeros((8, 12))
    times = np.zeros((8, 12))
    for i, j in product(range(0, 8), range(0, 12)):
        fig = plt.subplot(8, 12, i*12+j+1)
        data = plates_stack[:, i, j]
        if grad:
            current_speed = (float(d_time/np.sort(data)[-4]*60))
            # current_time = (d_time*np.argmax(data)-(d_time/np.max(data)*3))
            current_time = d_time*np.argsort(data)[-4]
            plt.title('%s, %s' % ('{0:.0f}'.format(current_speed),
                                  '{0:.2f}'.format(current_time)))
            speeds[i, j] = current_speed
            times[i, j] = current_time
        plt.plot(data.tolist())
        plt.ylim((np.min(plates_stack), np.max(plates_stack)))
        if j != 0:
            fig.set_yticklabels([])
        if i != 7:
            fig.set_xticklabels([])
        if i == 7:
            tick_lbls = ['{0:.1f}'.format(d_time*int(item)) for item in fig.get_xticks()]
            # tick_lbls = [(lambda x, term: term if x % 2 == 0 else '')(_i, val) for _i, val in enumerate(tick_lbls)]
            fig.set_xticklabels(tick_lbls)
            for tick in fig.xaxis.get_major_ticks():
                tick.label.set_fontsize(10)
                tick.label.set_rotation('vertical')

    if dumpType in ['raw', 'log', 'grad']:
        pickle.dump(plates_stack, open(Locations['dump'][dumpType], 'w'))
        plt.savefig(Locations['output'][dumpType], dpi=300)
        if dumpType == 'grad':
            pickle.dump((speeds, times), open(Locations['dump']['times'], 'w'))
    if NoShow:
        plt.clf()
    else:
        plt.show()


def group_plot(plates_stack, zoomlist):
    timepad = np.linspace(0, d_time*plates_stack.shape[0], num=plates_stack.shape[0])
    for sublist in zoomlist:
        legend = []
        for elt in sublist:
            plt.plot(timepad, plates_stack[:, elt[0], elt[1]])
            legend.append(str('%s %s') % (ascii_lowercase[elt[0]], elt[1] + 1))
        plt.legend(legend, loc='upper left')
        plt.savefig(os.path.join( base_folder_name, 'outputs/group_im_%s.png' % int(time() - tinit)))
        plt.show()
        # plt.clf()


def analyse(plates_stack, zoomlist, NoShow=False):
    if intermediate_show:
        plot_growth(plates_stack, dumpType='raw', NoShow=NoShow)
    reference_std = np.std(plates_stack[:, 0, 0])*2
    print reference_std
    log_stack = np.log10(plates_stack)/np.log10(2)
    for i, j in product(range(0, 8), range(0, 12)):
        log_stack[:, i, j] = log_stack[:, i, j] - np.mean(log_stack[range(0, 3), i, j])
    if intermediate_show:
        plot_growth(log_stack, dumpType='log', NoShow=NoShow)
    grad_stack = np.zeros(log_stack.shape)
    for i, j in product(range(0, 8), range(0, 12)):
        grad_stack[:, i, j] = np.gradient(gf_1d(log_stack[:, i, j], 2))
    if intermediate_show:
        plot_growth(grad_stack, True, dumpType='grad', NoShow=NoShow)
    group_plot(plates_stack, zoomlist)
    group_plot(log_stack, zoomlist)
    group_plot(grad_stack, zoomlist)


def correct(plate, position, injections):
    new_plate = np.zeros((plate.shape[0]+injections-1, plate.shape[1], plate.shape[2]))
    new_plate[:position+1, :, :] = plate[:position+1, :, :]
    new_plate[position+injections-1:, :, :] = plate[position:, :, :]
    diffplate = (plate[position+1, :, :] - plate[position, :, :]) / float(injections)
    for i in range(1, injections):
        new_plate[position+i, :, :] = plate[position, :, :] + diffplate * i
    return new_plate


def gaussian_process_regress(timeseries, std, timestamps=None, show=False):

    def show_routine():
        plt.figure()

        ax1 = plt.subplot(2, 2, 1)
        plt.title('OD measurement')
        plt.errorbar(timestamps.ravel(), timeseries, errors, fmt='r.', markersize=10, label=u'Observations')
        plt.plot(pre_timestamps, y_pred, 'b-', label=u'Regression')
        plt.fill(np.concatenate([pre_timestamps, pre_timestamps[::-1]]),
                                np.concatenate([y_pred - 1.9600 * sigma,
                                (y_pred + 1.9600 * sigma)[::-1]]),
                                alpha=.5, fc='b', ec='None', label='95% confidence interval')
        plt.xlabel('Time (hours)')
        plt.ylabel('OD')
        plt.legend(loc='upper left')

        plt.subplot(2, 2, 2, sharex=ax1)
        plt.title('Division speed')
        l2_growth = np.log(timeseries.ravel()[1:]/timeseries.ravel()[:-1])/np.log(2)/(timestamps[1:]-timestamps[:-1]).ravel()
        pred_l2_growth = np.log(y_pred.ravel()[1:]/y_pred.ravel()[:-1])/np.log(2)/(pre_timestamps[1:]-pre_timestamps[:-1]).ravel()
        plt.plot(timestamps[1:], l2_growth, 'r.',  label=u'Observations')
        plt.plot(pre_timestamps[1:], pred_l2_growth,  label=u'Regression')
        plt.yscale('log')
        plt.xlabel('Time (hours)')
        plt.ylabel('time to divide (minutes)')

        plt.subplot(2, 2, 3, sharex=ax1)
        plt.title('Instantaneous doubling time')
        l2_growth = 60/(np.log(timeseries.ravel()[1:]/timeseries.ravel()[:-1])/np.log(2)/(timestamps[1:]-timestamps[:-1]).ravel())
        pred_l2_growth = 60/(np.log(y_pred.ravel()[1:]/y_pred.ravel()[:-1])/np.log(2)/(pre_timestamps[1:]-pre_timestamps[:-1]).ravel())
        l2_growth[np.abs(l2_growth) > 300] = np.nan
        pred_l2_growth[np.abs(pred_l2_growth) > 300] = np.nan
        plt.plot(timestamps[1:], l2_growth, 'r.',  label=u'Observations')
        plt.plot(pre_timestamps[1:], pred_l2_growth,  label=u'Regression')
        plt.xlabel('Time (hours)')
        plt.ylabel('time to divide (minutes)')

        plt.subplot(2, 2, 4, sharex=ax1)
        plt.title('Cummulative doubling time')
        l2_growth = 60/(np.log(timeseries.ravel()[1:]/timeseries.ravel()[1])/np.log(2)/(timestamps[1:]-timestamps[1]).ravel())
        pred_l2_growth = 60/(np.log(y_pred.ravel()[1:]/y_pred.ravel()[1])/np.log(2)/(pre_timestamps[1:]-pre_timestamps[1]).ravel())
        l2_growth[np.abs(l2_growth) > 300] = np.nan
        pred_l2_growth[np.abs(pred_l2_growth) > 300] = np.nan
        plt.plot(timestamps[1:], l2_growth, 'r.',  label=u'Observations')
        plt.plot(pre_timestamps[1:], pred_l2_growth,  label=u'Prediction')
        plt.xlabel('$Time (hours)$')
        plt.ylabel('$division per hour$')

        plt.show()

    if timestamps is None:
        timestamps = np.linspace(0, timeseries.shape[0]*d_time, timeseries.shape[0])[:, np.newaxis]

    pre_timestamps = timestamps.copy()

    keep_mask = timeseries > 0.001
    timestamps = timestamps[:, 0][keep_mask][:, np.newaxis]
    timeseries = timeseries[keep_mask]

    nugget = np.convolve(timeseries, np.ones((5,))/5, mode='valid')
    nugget = np.lib.pad(nugget, (2, 2), 'edge')
    errors = np.abs(nugget - timeseries)
    errors[errors < std] = std
    nugget = np.power((errors), 2)

    gp = GaussianProcess(regr='linear', corr='squared_exponential', theta0=1,
                            thetaL=1e-2, thetaU=10,
                            nugget=nugget,
                            random_start=100)

    gp.fit(timestamps, timeseries)
    y_pred, MSE = gp.predict(pre_timestamps, eval_MSE=True)
    sigma = np.sqrt(MSE)

    if show:
        show_routine()

    elif np.any(y_pred < 0.001):
        # show_routine()
        pass

    return y_pred, sigma


def gaussian_process_wrapper(bulk_arguments):
    i, j, pl, std = bulk_arguments
    print 'loess', i, j
    if not debug:
        return (i, j), gaussian_process_regress(pl, std)
    if debug:
        if i < 1 or j < 1:
            return (i, j), gaussian_process_regress(pl, std)
        else:
            return (i, j), gaussian_process_regress(pl, std, show=True)


def map_adapter(plate, std):
     for i, j in product(range(0, 8), range(0, 12)):
         yield i, j, plate[:, i, j], std


def generate_reference_mask(plate_2):

    def extracted_growth_detector(_1D_array):
        return np.percentile(_1D_array, 99.5) - np.percentile(_1D_array, 0.5) < 2

    ref_mask = np.zeros((8, 12)).astype(np.bool)
    ref_mask[:, 0] = True
    ref_mask[:, 11] = True
    ref_mask[0, :] = True
    ref_mask[7, :] = True

    growth_not_detected = np.apply_along_axis(extracted_growth_detector, 0, plate_2)
    ref_mask = np.logical_and(ref_mask, growth_not_detected)
    medians = np.apply_along_axis(np.median, 0, plate_2)
    main_median = np.median(medians[ref_mask])
    means_not_off = np.abs(medians - main_median) < 0.2
    ref_mask = np.logical_and(ref_mask, means_not_off)
    timed_ref_mask = np.repeat(ref_mask[np.newaxis, :, :], plate_2.shape[0], axis=0)

    return timed_ref_mask


def fine_tune(plate_1, tune):
    # print 11, np.min(plate_1)
    plate_1[plate_1 < tune] = tune
    # print 12, np.min(plate_1)
    return plate_1


def loess(plate):
    refsample = plate_3D_array[generate_reference_mask(plate)]
    std = 10*np.std(refsample)

    # f_tune = np.percentile(refsample[refsample > 0], 1)
    # print 1, f_tune
    plate = plate - np.min(refsample)
    plate = fine_tune(plate, 0.0001)

    re_plate = plate.copy()

    # print 3,  np.min(re_plate)

    retset = map(gaussian_process_wrapper, map_adapter(plate, std))
    for ((i, j), (ret, _)) in retset:
        re_plate[:, i, j] = ret
    fine_tune(re_plate, 0.0001)

    # print 4,  np.min(re_plate)

    return re_plate


def smooth_and_interpolate(plate):
    refpoints = loess(plate)
    # deviation = np.abs(refpoints - plate)
    # _99 = np.percentile(deviation, 99)
    # plate[deviation > _99] = -1
    # refpoints = loess(plate)
    return refpoints


def del_exception(plate, position):
    """
    deletes a position and replaces it by an interpolation
    :param plate:
    :param position:
    :return:
    """
    new_plate = np.zeros((plate.shape[0]-1, plate.shape[1], plate.shape[2]))
    new_plate[:position, :, :] = plate[:position, :, :]
    print position
    new_plate[position:, :, :] = plate[position+1:, :, :]
    return new_plate


def del_range(plate, positionList):
    for position in sorted(positionList, reverse=True):
        plate = del_exception(plate, position)
    return plate


def correct_ODs(plate_3D_array):
    corrector = load_corrector()
    vect_corrector = np.vectorize(corrector)
    return  vect_corrector(plate_3D_array)


if __name__ == "__main__":
    path = os.path.join(file_location, file_name)
    pad_pth = os.path.join(pad_location, pad_name)
    intermediate_show = True
    read_pad(pad_pth)
    plate_3D_array = extract_plate_dit(path)
    plate_3D_array = del_exception(plate_3D_array, 18)
    plate_3D_array = del_exception(plate_3D_array, 11)
    plate_3D_array = correct_ODs(plate_3D_array)
    plate_3D_array = smooth_and_interpolate(plate_3D_array)

    # plate_3D_array = smooth_plate(plate_3D_array, 2)
    # plate_3D_array = del_exception(plate_3D_array, 220)
    # plate_3D_array = del_exception(plate_3D_array, 220)
    # plate_3D_array = correct(plate_3D_array, 219, 6)
    # del_range(plate_3D_array, range(220,222))
    zlist = []
    zlist = [
                # [(1, 1), (5, 1), (1, 3), (4, 3), (5, 3) ],
                # [(1, 1), (1, 3), (1, 6), (1, 10)],
                # [(3, 1), (3, 3), (3, 6), (3, 10)],
                # [(5, 1), (4, 3), (5, 3)],
                # [(2, 6), (6, 6), (3, 7), (4, 7), (1, 2), (2, 2)],
                ]
    analyse(plate_3D_array, zlist, NoShow=True)
