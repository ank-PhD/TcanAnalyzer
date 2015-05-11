__author__ = 'ank'

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from functools import partial
from adapters import hung_ji_adapter
from common import Plate
from configs import Locations
from csv import writer


def growth(timepoints, stepness, maxVal, midpoint, delay):
    # maxVal = 50  # works well
    # we will just remember all the parameters; this functions seems to be working better for some reason then the other ones.
    # if delay > 0: # delay control
    preliminary = maxVal/(1. + np.exp(-np.log(2)*stepness*(timepoints - midpoint)))
    # preliminary[timepoints < delay] = 0
    return preliminary


def minimize_function_builder2(timepoints, array):

    def minfunct(paramset):
        stepness, maxVal, midpoint, delay = paramset
        estimate = growth(timepoints, stepness, maxVal, midpoint, delay)
        difference = np.abs(estimate - array)
        return np.sum(difference**2)

    return minfunct


def flatten_and_group2(time_array, plate_3D):

    def iterative_fit(time, value_set):

        def gu(min, max):
            return np.random.uniform(min, max)

        def show(fit_params):
            higher_time = np.linspace(np.min(time), np.max(time), 100)
            plt.plot(time, value_set, 'r*')
            plt.plot(higher_time, growth(higher_time, 1/fit_params[0], *fit_params[1:-1]), 'k', label=' doubling time: %.2f h\n max: %.2f \n midpoint: %.0f h\n lag: %.0f h\n error: %.4f\n '% tuple(fit_params))
            plt.legend(loc='upper left', prop={'size':10})
            plt.show()

        v_set = np.array(value_set)
        v_set -= np.min(v_set)

        bounds = [(0.001, 0.99), (0.002, 0.1), (0, 120), (0, 1)] #TODO: growth-wise lag optimisaiton is skrewed

        ffit, errcode = fit_with_flat(time, v_set, bounds=bounds)
        if ffit[-1] > 0.65 and errcode != 1:
            for i in range(0, 10):
                start = [gu(*bound) for bound in bounds]
                ffit, errcode = fit_with_flat(time, v_set, start_point=start, bounds=bounds)
                if ffit[-1] < 0.4:
                    break
            if ffit[-1] > 0.65:
                errcode = 2

        print errcode
        if errcode == 2:
            print "errparams %s, %s" %(ffit, errcode)
            show(ffit)

        ffit.append(errcode)
        return ffit

    def fit_with_flat(time, v_set, start_point=[0.05, 0.005, 40., 5.], bounds=[(0.01, 0.99), (0.001, 0.5), (0, 100), (0, 10)] ):

        take_off = np.max(v_set[1:-1])
        if take_off < 0.005:
            return ['inf', 'NA', 'NA', 'NA', np.mean(np.abs(np.mean(v_set, axis=0)))], 1

        mfunct = minimize_function_builder2(np.array(time), v_set)
        OR_object = minimize(mfunct, start_point, bounds=bounds)
        popt = OR_object.x
        if OR_object.success:
            return [1./popt[0]] + popt[1:].tolist() + [np.mean(np.abs(growth(np.array(time), *popt)-v_set))/np.mean(np.abs(v_set))], 0

        else:
            print OR_object.message,
            print OR_object.x
            popt = OR_object.x
            return [1./popt[0]] + popt[1:].tolist() + [np.mean(np.abs(growth(np.array(time), *popt)-v_set))/np.mean(np.abs(v_set))], -1

    _1D_fit = partial(iterative_fit, time_array)
    result = np.apply_along_axis(_1D_fit, axis=0, arr=plate_3D)
    print result.shape
    print 'result: \n', result
    return Plate(result).export_values()


def write_resuts(_2D_list):
    with open(Locations['output']['hjt'], 'wb') as destination_file:
        destination = writer(destination_file)
        destination.writerows(_2D_list)


if __name__ == "__main__":
    hjt_dict = hung_ji_adapter()
    accumulator = [['row', 'column', 'hours to double', 'maximum value', 'midpoint', 'delay (to be ignored)']]
    for value, (time_, array_3D_) in hjt_dict.iteritems():
        accumulator.append(['Concentration', value])
        accumulator += flatten_and_group2(time_, array_3D_)
    write_resuts(accumulator)
