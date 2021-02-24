DMSP Dispersion Event Identification Software
---------------------------------------------

This repository holds code and notebooks associated with the automated dispersion event search for use with data from the DMSP satellites. The main scripts are:

* `download.sh` - Download data using [MagridalWeb](http://cedar.openmadrigal.org/madrigalDownload). Edit the `yyyy` variable in the header to specify the year, and adjust the dates as needed.
* `search_dispersion_events.py` - Script to search for events-- you pass it a DMSP hdf5 filename and it prints a table of discovered events. You can also call functions in this module directly to get the list of events back as a Pandas DataFrame.


This code was developed by Daniel da Silva at NASA Goddard Spaceflight Center, who may be contacted at [mail@danieldasilva.org](mailto:mail@danieldasilva.org), [daniel.e.dasilva@nasa.gov](mailto:daniel.e.dasilva@nasa.gov), or [ddasilva@ursa.edu](mailto:ddasilva@usra.edu).

