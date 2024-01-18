#!/usr/bin/env python
"""Writes a case file for use with analyzing a storm.

See also: run_model.py
"""
import argparse
import glob
import json
import os
from termcolor import cprint


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('run_name')
    parser.add_argument('sc_num', type=int, help="Spacecraft Number (e.g., 16)")
    parser.add_argument('--reverse-effect', action='store_true')
    parser.add_argument('--inverse-effect', action='store_true')
    args = parser.parse_args()
    
    # Edit here
    reverse_effect = args.reverse_effect
    inverse_effect = args.inverse_effect
    
    plot_output = f'output/{args.run_name}_F{args.sc_num}/plots'
    event_output = f'output/{args.run_name}_F{args.sc_num}/{args.run_name}_F{args.sc_num}.csv'

    # Get OMNIweb files
    omniweb_glob = 'data/' + args.run_name + '/**/omni*.cdf'
    omniweb_files = []
    omniweb_files.extend(glob.glob(omniweb_glob, recursive=True))

    # Get DMSP files
    dmsp_flux_glob = 'data/' + args.run_name + f'/Satellite_F{args.sc_num}/**/*e.*.hdf5'
    dmsp_flux_files = []
    dmsp_flux_files.extend(glob.glob(dmsp_flux_glob, recursive=True))
    dmsp_flux_files.sort()
    
    dmsp_magn_glob = 'data/' + args.run_name + f'/Satellite_F{args.sc_num}/**/*s1.*.hdf5'
    dmsp_magn_files = []
    dmsp_magn_files.extend(glob.glob(dmsp_magn_glob, recursive=True))
    dmsp_magn_files.sort()

    # Make plot output dir if does not exist
    os.makedirs(plot_output, exist_ok=True)
    
    # Write case fiel
    case_file = {
        'STORM_NAME': args.run_name,
        'DMSP_FLUX_FILES': dmsp_flux_files,
        'DMSP_MAGN_FILES': dmsp_magn_files,
        'OMNIWEB_FILES': omniweb_files,
        'PLOT_OUTPUT': plot_output,
        'EVENT_OUTPUT': event_output,
        'REVERSE_EFFECT': reverse_effect,
        'INVERSE_EFFECT': inverse_effect,
    }

    case_filename = f'case_files/{args.run_name}_F{args.sc_num}.json'
    fh = open(case_filename, 'w')
    json.dump(case_file, fh, indent=4)
    fh.write('\n')
    fh.close()
    
    cprint(f'Wrote case file to path {case_filename}', 'green')
    cprint(f"Plots will be written to {plot_output}", "yellow")
    cprint(f"CSV output will be written to {event_output}", "yellow")


if __name__ == '__main__':
    main()
