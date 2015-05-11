__author__ = 'ank'

import numpy as np
from pprint import pprint
from itertools import product

class Plate(object):

    def __init__(self, _2D_plate=None):
        if _2D_plate is None:
            self.np_array = np.zeros((8, 12))
        else:
            self.np_array = _2D_plate

        self.r_map = lambda x: ord(x) - 65
        self.c_map = lambda x: int(x) - 1

        self.r_rev = lambda x: chr(x + 65)
        self.c_rev = lambda x: int(x) + 1

    def set_elt(self, row, column, value):
        self.np_array[self.r_map(row), self.c_map(column)] = float(value)

    def export_values(self):
        accumulator = []
        for i, j in product(range(0, 8), range(0, 12)):
            t_list = self.np_array[:, i, j].tolist()
            t_list.insert(0, self.c_rev(j))
            t_list.insert(0, self.r_rev(i))
            accumulator.append(t_list)
        return accumulator


if __name__ == "__main__":
    plt = Plate()
    plt.set_elt('A', 1, 4.0)
    print plt.np_array
    pprint(plt.export_values())