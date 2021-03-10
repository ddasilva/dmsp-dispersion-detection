#!/usr/bin/env python
"""Writes a case file for use with analyzing a storm.

See also: visualize_storm.py
"""

import glob
import json
import os
from termcolor import cprint


def main():
    # Edit here
    storm_name = 'Nov20_2003_storm'
    omniweb_file = 'data/' + storm_name + '/omni_hro_1min_20031101_v01.cdf'
    plot_output = 'plots/' + storm_name
    event_output = 'data/' + storm_name + '.csv'

    days = list(range(18, 24))
    dmsp_glob = 'data/'+ storm_name + '/**/*200311%02d*.hdf5'

    # Get DMSP files
    dmsp_files = []
    for day in days:
        dmsp_files.extend(glob.glob(dmsp_glob%day, recursive=True))

    dmsp_files.sort()
    
    # Make plot output dir if does not exist
    os.makedirs(plot_output, exist_ok=True)
    
    # Write case fiel
    case_file = {
        'STORM_NAME': storm_name,
        'DMSP_FILES': dmsp_files,
        'OMNIWEB_FILE': omniweb_file,
        'PLOT_OUTPUT': plot_output,
        'EVENT_OUTPUT': event_output,
        
    }
    case_filename = 'case_files/' + storm_name + '.json'
    fh = open(case_filename, 'w')
    json.dump(case_file, fh, indent=4)
    fh.write('\n')
    fh.close()

    cprint('Wrote case file to path ' + case_filename, 'green')

if __name__ == '__main__':
    main()
