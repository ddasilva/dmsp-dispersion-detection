"""Helper script to download DMSP data from MadrigalWeb"""

from dateutil.relativedelta import relativedelta
from datetime import datetime, timedelta
import subprocess
import argparse
import os
from termcolor import cprint
import urllib.request

parser = argparse.ArgumentParser()
parser.add_argument('start_time', help='Start time, in mm/dd/yyyy')
parser.add_argument('end_time', help='End time, in mm/dd/yyyy')
parser.add_argument('run_name', help='Name of run. Data gets put in data/$run_name')
args = parser.parse_args()

toks = args.start_time.split('/')
start_date = datetime(int(toks[2]), int(toks[0]), int(toks[1])) - timedelta(days=1)
toks = args.end_time.split('/')
end_date = datetime(int(toks[2]), int(toks[0]), int(toks[1]))

cur_date = start_date

while cur_date < end_date:
    month, year = cur_date.month, cur_date.year
    cprint(f"Downloading omniweb data for {month}/{year} from URL:", "green")    

    fname = f"omni_hro_1min_{year}{month:02d}01_v01.cdf"
    url = f"https://spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/{year}/{fname} "
    outfile = f"data/{args.run_name}/omni/{year}/{fname}"

    cprint(f"  {url}", "cyan")
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    urllib.request.urlretrieve(url, outfile)
    assert os.path.exists(outfile)    
    cprint(f"  ...and writing to {outfile}", "white")
    #subprocess.check_call(cmd, shell=True)
    cur_date += relativedelta(months=1)
