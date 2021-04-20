#!/bin/bash

# Download DMSP data from Madrigal
for i in $(seq 16 19); do
    globalDownload.py --verbose --url=http://cedar.openmadrigal.org --outputDir=./data/Dec20_2015_storm/Satellite_F${i} --user_fullname='Daniel+da+Silva' --user_email=daniel.e.dasilva@nasa.gov --user_affiliation='NASA' --format='hdf5' --startDate="12/18/2015" --endDate="12/23/2015" --inst=8100 --kindat=102${i},101${i} &
done

wait  # wait for all background processes to complete

# Download IMF data from OMNIWeb
wget https://spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/2015/omni_hro_1min_20151201_v01.cdf -O data/Dec20_2015_storm/omni_hro_1min_20151201_v01.cdf


