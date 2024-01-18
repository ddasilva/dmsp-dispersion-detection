#!/bin/bash

# Download data from Madrigal. Requires the Python module to be installed,
# which can be obtained from http://cedar.openmadrigal.org/madrigalDownload

# Edit to specify year. Edit the --start-date and --end-date below
# to further adjust the range.

if [ "$#" -ne 3 ]; then
    echo "Usage: sh $0 <mm/dd/yyyy> <mm/dd/yyyy> <run_name> # start time and end time"
    exit
fi

start_time=$1
end_time=$2
run_name=$3

echo start_time=$start_time
echo end_time=$end_time
echo run_name=$run_name


for i in $(seq 16 19); do   
    globalDownload.py --verbose --url=http://cedar.openmadrigal.org --outputDir=./data/$run_name/Satellite_F${i} --user_fullname=$(whoami) --user_email=daniel.e.dasilva@nasa.gov --user_affiliation='NASA' --format='hdf5' --startDate="$start_time" --endDate="$end_time" --inst=8100 --kindat=102${i} &
done

wait
