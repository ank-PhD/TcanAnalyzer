import os
import numpy as np
from csv import reader
from collections import defaultdict
from common import Plate
from pprint import pprint


def hung_ji_adapter():
    folder_ = 'C:/Users/ank/Desktop'
    file_ = 'HJT_fittness.csv'
    dose_curves = {}
    with open(os.path.join(folder_, file_)) as source_file:
        source_ = reader(source_file)
        source_.next()
        dose2data_frame = defaultdict(lambda:defaultdict(Plate))
        for line in source_:
            dose, time, col, row, intensity = line
            plate = dose2data_frame[dose][time]
            plate.set_elt(row, col, intensity)
        for dose, timedict in dose2data_frame.iteritems():
            bound_tuple = [[], []]
            dose_curves[dose] = bound_tuple
            for time in sorted(timedict.keys()):
                bound_tuple[0].append(time)
                bound_tuple[1].append(timedict[time].np_array)
            bound_tuple[0] = np.array(bound_tuple[0]).astype(np.float)
            bound_tuple[1] = np.array(bound_tuple[1]).astype(np.float)
    return dose_curves


if __name__ == "__main__":
    pprint(hung_ji_adapter())