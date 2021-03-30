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
    storm_name = 'Long_Term_Trend'
    plot_output = 'plots/' + storm_name
    event_output = 'data/' + storm_name + '.csv'
    reverse_effect = False
    
    # Get OMNIweb files
    omniweb_glob = 'data/' + storm_name + '/omniweb/**/*.cdf'
    omniweb_files = []
    omniweb_files.extend(glob.glob(omniweb_glob, recursive=True))

    # Get DMSP files
    dmsp_glob = 'data/' + storm_name + '/Satellite_*/**/*.cdf'
    dmsp_files = []
    dmsp_files.extend(glob.glob(dmsp_glob, recursive=True))

    dmsp_files.sort()
    
    # Make plot output dir if does not exist
    os.makedirs(plot_output, exist_ok=True)
    
    # Write case fiel
    case_file = {
        'STORM_NAME': storm_name,
        'DMSP_FILES': dmsp_files,
        'OMNIWEB_FILES': omniweb_files,
        'PLOT_OUTPUT': plot_output,
        'EVENT_OUTPUT': event_output,
        'REVERSE_EFFECT': reverse_effect,
    }

    case_filename = 'case_files/' + storm_name + '.json'
    fh = open(case_filename, 'w')
    json.dump(case_file, fh, indent=4)
    fh.write('\n')
    fh.close()

    cprint('Wrote case file to path ' + case_filename, 'green')

if __name__ == '__main__':
    main()
