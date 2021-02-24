#!/bin/bash

# Download data from Madrigal. Requires the Python module to be installed,
# which can be obtained from http://cedar.openmadrigal.org/madrigalDownload

# Edit to specify year. Edit the --start-date and --end-date below
# to further adjust the range.
yyyy=2010

for i in $(seq 16 19); do
    
    /home/daniel/anaconda3/bin/globalDownload.py --verbose --url=http://cedar.openmadrigal.org --outputDir=./data/validate/$yyyy/Satellite_F${i} --user_fullname='Daniel+da+Silva' --user_email=daniel.e.dasilva@nasa.gov --user_affiliation='NASA' --format='hdf5' --startDate="1/1/$yyyy" --endDate="12/31/$yyyy" --inst=8100 --kindat=102${i} &
done

wait
