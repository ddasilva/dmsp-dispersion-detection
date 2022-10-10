DMSP Dispersion Event Identification Software
---------------------------------------------

This repository holds code and notebooks associated with the automated dispersion event search for use with data from the DMSP satellites. The main scripts are:

* `download_dmsp.sh` - Download data using [MagridalWeb](http://cedar.openmadrigal.org/madrigalDownload). Edit the `yyyy` variable in the header to specify the year, and adjust the dates as needed.
* `download_omniweb.sh` - Downoad data from OMNIWeb HTTP Server. Edit the `startyear` and `endyear` variables as needed.
* `run_model.py` - Script to search data and plot discovered events-- you pass it a case file created with `make_case_file.py`.
* `make_case_file.py` - Generate a case file which holds all the inputs for `run_model.py`. Edit the variables at the top of the code and re-run the script to get it to write a file.

Requirements
* [Python 3 with Miniconda/Anaconda](https://docs.conda.io/en/latest/miniconda.html) - programming language, see `environment.yml` for module dependencies.
* [Wget](https://www.gnu.org/software/wget/) - recent version needed, required for downloading data from web and FTP.

This code was developed by Daniel da Silva at NASA Goddard Spaceflight Center, who may be contacted at [mail@danieldasilva.org](mailto:mail@danieldasilva.org), [daniel.e.dasilva@nasa.gov](mailto:daniel.e.dasilva@nasa.gov), or [ddasilva@ursa.edu](mailto:ddasilva@usra.edu).

