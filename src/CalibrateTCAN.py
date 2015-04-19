__author__ = 'ank'

import numpy as np
from matplotlib import pyplot as plt
from scipy import odr
from csv import reader

errors_location = '../normalize/Errors.csv'
alignement_location = '../normalize/Values.csv'

def lin(p, x):
    a, b, c, d = p
    return c*x

def pol2(p, x):
     a, b, c, d = p
     return b*x**2 + c*x + d


def pol3(p, x):
    a, b, c, d = p
    return a*x**3 + b**x**2 + c*x + d


def regress(x, y, x_sd, y_sd, function_to_fit=pol3, name_to_plot='', figure_no=1):

    def plot_result():
        x_fit = np.linspace(np.min(x)*0.95, np.max(x)*1.05, 1000)
        y_fit = function_to_fit(out.beta, x_fit)
        lin_fit = lin(lin_out.beta, x_fit)
        plt.subplot(2, 2, figure_no)
        plt.title(name_to_plot+': \n %.2fx^3 + %.2fx^2 + %.2fx + %.2f v.s. %.2fx. \n Res var gain: x %.2f' % tuple(out.beta.tolist()+[lin_out.beta.tolist()[2]]+[lin_out.res_var/out.res_var]))
        plt.errorbar(x, y, xerr=x_sd, yerr=y_sd, linestyle='None', marker='x')
        plt.plot(x_fit, y_fit, 'g')
        plt.plot(x_fit, lin_fit, 'r')
        plt.autoscale(tight=True)

    model = odr.Model(function_to_fit)
    data = odr.RealData(x, y, sx=x_sd, sy=y_sd)
    _odr = odr.ODR(data, model, beta0=[1., 1., 10., 0.01])
    out = _odr.run()

    lin_model = odr.Model(lin)
    lin_odr = odr.ODR(data, lin_model, beta0=[0., 0., 10., 0.01])
    lin_out = lin_odr.run()

    lin_out.pprint()

    plot_result()

    return out.beta


def read_normalization_tables(limit = 15):

    def normalize_column(array_1D):
        return array_1D - min(array_1D)

    accumulator = []
    names = []
    with open(errors_location, 'rb') as csvfile:
        rdr = reader(csvfile)
        names = rdr.next()
        for row in rdr:
            accumulator.append(row)
    errors = np.array(accumulator).astype(np.float)

    accumulator = []
    with open(alignement_location, 'rb') as csvfile:
        rdr = reader(csvfile)
        rdr.next()
        for row in rdr:
            accumulator.append(row)
    correspondence_tables = np.array(accumulator).astype(np.float)

    correspondence_tables = np.apply_along_axis(normalize_column, axis=0, arr=correspondence_tables)

    msk = correspondence_tables[:, 0] < limit
    errors = errors[msk, :]
    correspondence_tables = correspondence_tables[msk, :]

    print correspondence_tables.shape

    return names, errors.T, correspondence_tables.T


def total_regression(errors, correspondance_tables, names):
    ref_err, ref_reads = errors[0, 0], correspondance_tables[0, :]
    for _i in range(1, errors.shape[0]):
        cor_err, cor_reads  = errors[_i, 0], correspondance_tables[_i, :]
        regress(cor_reads, ref_reads, cor_err, ref_err, name_to_plot=names[_i], figure_no=_i)


if __name__ == "__main__":
    names, err, table = read_normalization_tables(3)
    total_regression(err, table, names)
    plt.show()
