#!/bin/bash

# Download DMSP data from Madrigal
for i in $(seq -w 13 16); do
    globalDownload.py --verbose --url=http://cedar.openmadrigal.org --outputDir=./data/Nov20_2003_storm/Satellite_F${i} --user_fullname='Daniel+da+Silva' --user_email=daniel.e.dasilva@nasa.gov --user_affiliation='NASA' --format='hdf5' --startDate="11/18/2003" --endDate="11/23/2003" --inst=8100 --kindat=102${i} &
done

wait  # wait for all background processes to complete

# Download IMF data from OMNIWeb
wget https://spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/2003/omni_hro_1min_20031101_v01.cdf -O data/Nov20_2003_storm/omni_hro_1min_20031101_v01.cdf


