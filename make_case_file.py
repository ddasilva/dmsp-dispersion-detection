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
    storm_name = 'LJ2015_Sep_to_Dec'
    sat_num = 19
    reverse_effect = False
    inverse_effect = True
    
    plot_output = 'plots/' + storm_name
    event_output = f'data/{storm_name}_F{sat_num}.csv'
    
    # Get OMNIweb files
    omniweb_glob = 'data/' + storm_name + '/**/omni*.cdf'
    omniweb_files = []
    omniweb_files.extend(glob.glob(omniweb_glob, recursive=True))

    # Get DMSP files
    dmsp_flux_glob = 'data/' + storm_name + f'/Satellite_F{sat_num}/**/*e.*.hdf5'
    dmsp_flux_files = []
    dmsp_flux_files.extend(glob.glob(dmsp_flux_glob, recursive=True))
    dmsp_flux_files.sort()
    
    dmsp_magn_glob = 'data/' + storm_name + f'/Satellite_F{sat_num}/**/*s1.*.hdf5'
    dmsp_magn_files = []
    dmsp_magn_files.extend(glob.glob(dmsp_magn_glob, recursive=True))
    dmsp_magn_files.sort()

    # Make plot output dir if does not exist
    os.makedirs(plot_output, exist_ok=True)
    
    # Write case fiel
    case_file = {
        'STORM_NAME': storm_name,
        'DMSP_FLUX_FILES': dmsp_flux_files,
        'DMSP_MAGN_FILES': dmsp_magn_files,
        'OMNIWEB_FILES': omniweb_files,
        'PLOT_OUTPUT': plot_output,
        'EVENT_OUTPUT': event_output,
        'REVERSE_EFFECT': reverse_effect,
        'INVERSE_EFFECT': inverse_effect,
    }

    case_filename = f'case_files/{storm_name}_F{sat_num}.json'
    fh = open(case_filename, 'w')
    json.dump(case_file, fh, indent=4)
    fh.write('\n')
    fh.close()

    cprint('Wrote case file to path ' + case_filename, 'green')


if __name__ == '__main__':
    main()
