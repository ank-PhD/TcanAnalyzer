__author__ = 'ank'

import os
base_folder_name = os.path.abspath(os.path.join(os.path.dirname(__file__), '.'))

Locations = {'dump':{
                    'raw': os.path.join(base_folder_name, 'dumps/raw.dmp'),
                    'log': os.path.join(base_folder_name, 'dumps/log.dmp'),
                    'grad': os.path.join(base_folder_name, 'dumps/grad.dmp'),
                    'times': os.path.join(base_folder_name, 'dumps/times.dmp')
                    },
             'output':{
                    'raw': os.path.join(base_folder_name, 'outputs/raw.png'),
                    'log': os.path.join(base_folder_name, 'outputs/log.png'),
                    'grad': os.path.join(base_folder_name, 'outputs/grad.png'),
                    'plate': os.path.join(base_folder_name, 'outputs/plate.png')
                    },
             'pad': os.path.join(base_folder_name, 'dumps/pad.dmp')
             }
