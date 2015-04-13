__author__ = 'ank'

import numpy as np
from matplotlib import pyplot as plt
from scipy import odr
from csv import reader

errors_location = '../normalize/Errors_alignement.csv'
alignement_location = '../normalize/Exponential_alignement.csv'

def pol2(p, x):
     a, b, c, d = p
     return b*x**2 + c*x + d


def pol3(p, x):
    a, b, c, d = p
    return a*x**3 + b**x**2 + c*x + d


def regress(x, y, x_sd, y_sd, function_to_fit=pol3, name_to_plot='', figure_no=1):

    def plot_result():
        x_fit = np.linspace(x[0], x[-1], 1000)
        y_fit = function_to_fit(out.beta, x_fit)
        plt.subplot(2, 2, figure_no)
        plt.title(name_to_plot+': %.2fx^3 + %.2fx^2 + %.2fx + %.2f' % tuple(out.beta.tolist()))
        plt.errorbar(x, y, xerr=x_sd, yerr=y_sd, linestyle='None', marker='x')
        plt.plot(x_fit, y_fit)

    model = odr.Model(function_to_fit)
    data = odr.RealData(x, y, sx=x_sd, sy=y_sd)
    _odr = odr.ODR(data, model, beta0=[0.0, 0.1, 1.0, 0.1])
    out = _odr.run()

    out.pprint()

    plot_result()

    return out.beta


def read_normalization_tables():
    accumulator = []
    names = []
    with open(errors_location, 'rb') as csvfile:
        rdr = reader(csvfile)
        names = rdr.next()
        accumulator.append(rdr.next())
    errors = np.array(accumulator).astype(np.float)

    accumulator = []
    with open(alignement_location, 'rb') as csvfile:
        rdr = reader(csvfile)
        rdr.next()
        for row in rdr:
            accumulator.append(row)
    correspondence_tables = np.array(accumulator).astype(np.float)

    return names, errors.T, correspondence_tables.T


def total_regression(errors, correspondance_tables, names):
    ref_err, ref_reads = errors[0, 0], correspondance_tables[0, :]
    for _i in range(1, errors.shape[0]):
        cor_err, cor_reads  = errors[_i, 0], correspondance_tables[_i, :]
        regress(cor_reads, ref_reads, cor_err, ref_err, name_to_plot=names[_i], figure_no=_i)


if __name__ == "__main__":
    names, err, table = read_normalization_tables()
    total_regression(err, table, names)
    plt.show()
