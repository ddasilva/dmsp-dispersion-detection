DMSP Dispersion Event Identification Software
---------------------------------------------

This repository holds code and notebooks associated with the automated dispersion event search for use with data from the DMSP satellites. The main scripts are:

* `download_dmsp.sh` - Download data using [MagridalWeb](http://cedar.openmadrigal.org/madrigalDownload). Edit the `yyyy` variable in the header to specify the year, and adjust the dates as needed.
* `download_omniweb.sh` - Downoad data from OMNIWeb HTTP Server. Edit the `startyear` and `endyear` variables as needed.
* `visualize_storm.py` - Script to search data and plot discovered events-- you pass it a case file created with `util_make_case_file.py`.
* `util_make_case_file.py` - Generate a case file which holds all the inputs for `visualize_storm.py`. Edit the variables at the top of the code and re-run the script to get it to write a file.

This code was developed by Daniel da Silva at NASA Goddard Spaceflight Center, who may be contacted at [mail@danieldasilva.org](mailto:mail@danieldasilva.org), [daniel.e.dasilva@nasa.gov](mailto:daniel.e.dasilva@nasa.gov), or [ddasilva@ursa.edu](mailto:ddasilva@usra.edu).

