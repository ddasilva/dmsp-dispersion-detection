#!/bin/bash


yyyy=2015

for i in $(seq 16 19); do
    globalDownload.py --verbose --url=http://cedar.openmadrigal.org --outputDir=./data/Dec20_2015_storm/Satellite_F${i} --user_fullname='Daniel+da+Silva' --user_email=daniel.e.dasilva@nasa.gov --user_affiliation='NASA' --format='hdf5' --startDate="12/1/$yyyy" --endDate="12/31/$yyyy" --inst=8100 --kindat=102${i} &
done

wait
