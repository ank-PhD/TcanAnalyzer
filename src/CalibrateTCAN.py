import numpy as np
from matplotlib import pyplot as plt
from scipy import odr
from csv import reader
import pickle
from configs import Locations


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
    re_switch = np.max(correspondence_tables[msk], axis=0)
    # print re_switch
    errors = errors[msk, :]
    correspondence_tables = correspondence_tables[msk, :]

    # print correspondence_tables.shape

    return names, errors.T, correspondence_tables.T, re_switch


def total_regression(errors, correspondance_tables, names):
    ref_err, ref_reads = errors[0, 0], correspondance_tables[0, :]
    param_table = []
    for _i in range(1, errors.shape[0]):
        cor_err, cor_reads  = errors[_i, 0], correspondance_tables[_i, :]
        param_table.append(regress(cor_reads, ref_reads, cor_err, ref_err, name_to_plot=names[_i], figure_no=_i))
    return param_table


def build_correction_function(corr_function, corr_params_a, corr_params_l, switch):

    def corrective_function(x):
        if x < switch:
            cf = corr_function(list(corr_params_l), x)
            if cf < 0.001:
                return 0.001
            else:
                return cf
        else:
            cf_l = corr_function(list(corr_params_l), x)
            cf_h = corr_function(list(corr_params_a), x)
            if cf_h < cf_l:
                return cf_l
            else:
                return cf_h

    # def corrective_function(x):
    #     cf = corr_function(list(corr_params_a), x)
    #     if cf < 0.005:
    #         return 0.005
    #     else:
    #         return cf

    return corrective_function


def create_corrector(switch=3, option=1):
    names_l, err_l, table_l, re_switch = read_normalization_tables(switch)
    reg_our_dynamic_l = total_regression(err_l, table_l, names_l)[option]
    plt.show()
    names, err, table, _ = read_normalization_tables()
    reg_our_dynamic = total_regression(err, table, names)[option]
    plt.show()
    pickle.dump((reg_our_dynamic_l, reg_our_dynamic, re_switch[option]), open(Locations['dump']['corrfunct'], 'w'))


def load_corrector(show=False):

    def show_correction():
        plt.title('OD correcting function')
        x = np.linspace(-0.01, 1.5, 500)
        plt.plot(x, np.vectorize(corrfunction)(x))
        plt.show()

    coeff_l, coeff_a, re_switch = pickle.load(open(Locations['dump']['corrfunct'], 'r'))

    corrfunction = build_correction_function(pol3, coeff_a, coeff_l, re_switch)

    if show:
        show_correction()

    return corrfunction


if __name__ == "__main__":
    # names, err, table = read_normalization_tables()
    # reg_our_dynamic = total_regression(err, table, names)[1]
    # plt.show()
    # build_correction_function(pol3, reg_our_dynamic)
    create_corrector(3, 1)
    load_corrector(True)


