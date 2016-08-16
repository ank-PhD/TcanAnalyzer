from chiffatools.high_level_os_methods import safe_dir_create
import os

base_folder_name = os.path.abspath(os.path.join(os.path.dirname(__file__), '.'))

dumps_folder = os.path.join(base_folder_name, 'dumps')
outputs_folder = os.path.join(base_folder_name, 'outputs')

safe_dir_create(dumps_folder)
safe_dir_create(outputs_folder)

Locations = {'dump':{
                    'raw': os.path.join(dumps_folder, 'raw.dmp'),
                    'log': os.path.join(dumps_folder, 'log.dmp'),
                    'grad': os.path.join(dumps_folder, 'grad.dmp'),
                    'times': os.path.join(dumps_folder, 'times.dmp'),
                    'corrfunct': os.path.join(dumps_folder, 'corrfunct.dmp')
                    },
             'output':{
                    'raw': os.path.join(outputs_folder, 'raw.png'),
                    'log': os.path.join(outputs_folder, 'log.png'),
                    'grad': os.path.join(outputs_folder, 'grad.png'),
                    'plate': os.path.join(outputs_folder, 'plate.png'),
                    'hjt': os.path.join(outputs_folder, 'hjt.csv')
                    },
             'pad': os.path.join(dumps_folder, 'pad.dmp')
             }
