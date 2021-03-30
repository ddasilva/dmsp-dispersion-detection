#!/bin/bash

# Download data Kp index data from GFZ Helmholtz Center Potsdam, which is
# located on the web at https://www.gfz-potsdam.de/en/kp-index/
#
# This downloads a text file which includes all Kp data since 1932.

outname="data/Kp_ap_since_1932.txt"

rm -f ${outname}

wget ftp://anonymous@ftp.gfz-potsdam.de/pub/home/obs/Kp_ap_Ap_SN_F107/Kp_ap_since_1932.txt -O ${outname}
