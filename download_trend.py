#!/bin/env python
"""Download DMSP data from the SPDF (Space Physics Data Facility) at GSFC.

This data is in CDF format, and is only as recent as 2014. To adjust the years 
and satellites, edit the header.
"""
import glob
import os
import subprocess

start_year = 2010
end_year = 2014
satellies = [16, 17, 18]


def download_omniweb():
    """Download OMNIWeb data from the SPDF"""
    # Transfer and oraganize data
    for year in range(start_year, end_year + 1):
        for mon in range(1, 13):
                outdir = f'data/Long_Term_Trend/omniweb/{year}/{mon:02d}'
                os.makedirs(outdir, exist_ok=True)
                
                fname = f'omni_hro_1min_{year}{mon:02d}01_v01.cdf'
                outpath = os.path.join(outdir, fname)
                cmd = f'wget -nc https://spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/{year}/{fname} -O {outpath}'
                print(f'Running cmd: {cmd}')
                os.system(cmd)


def download_dmsp():
    """Download DMSP data from the SPDF"""
    # Transfer Data
    for year in range(start_year, end_year + 1):
        for sat in satellies:
            cmd = (f'wget -r -np -nc -nd -l 1 -A cdf '
                   f'https://spdf.gsfc.nasa.gov/pub/data/dmsp/dmspf{sat}/ssj/precipitating-electrons-ions/{year}/')
            print(f'Running cmd: {cmd}')
            subprocess.check_output(cmd, shell=True)
            
    # Organize files 
    for fname in glob.glob('dmsp-*.cdf'):
        toks = fname.split('_')
        date = toks[-2]
        yyyy = date[:4]
        mm = date[4:6]
        sat = toks[0][-3:].upper()
    
        destdir = f'data/Long_Term_Trend/Satellite_{sat}/{yyyy}/{mm}/'
        dest = os.path.join(destdir, fname)    
        os.makedirs(destdir, exist_ok=True)
        
        print(f'Moving File {fname} to {dest}')
        os.rename(fname, dest)


def main():
    download_omniweb()
    download_dmsp()
    print('All done!')


if __name__ == '__main__':
    main()
