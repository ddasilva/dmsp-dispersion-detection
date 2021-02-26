#!/bin/bash

# This script downloads minutely magnetic field data from OMNIWeb between $startyear and $endyear,
# and saves the data only if the directory prior exists (this is to be run after the DMSP data is
# downloaded and creates the directories). This is so the data is only downloaded if there is DMSP
# data to match with it.

startyear=2011
endyear=2020

for year in $(seq ${startyear} ${endyear}); do
    for mon in $(seq -w 1 12); do	
	outdir=data/${year}/${mon}
	fname=omni_hro_1min_${year}${mon}01_v01.cdf
	[ -d ${outdir} ] && wget https://spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/$year/${fname} -O ${outdir}/${fname}
    done
done
