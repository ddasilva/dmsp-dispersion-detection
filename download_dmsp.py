"""Helper script to download DMSP data from MadrigalWeb"""

import subprocess
import argparse
import os
import shutil
from termcolor import cprint

parser = argparse.ArgumentParser()
parser.add_argument('start_time', help='Start time, in mm/dd/yyyy')
parser.add_argument('end_time', help='End time, in mm/dd/yyyy')
parser.add_argument('run_name', help='Name of run. Data gets put in data/$run_name')
parser.add_argument('--spacecraft_csv', default='16,17,18,19', help='CSV list of DMSP spacecraft numbers')
args = parser.parse_args()

sc_nums = args.spacecraft_csv.split(",")

if os.name == "nt":
    cmd_path = subprocess.getoutput("where globalDownload.py") 
else:
    cmd_path = subprocess.getoutput("which globalDownload.py")

cprint(f"Download data for spacecrafts: {sc_nums}", "green")

for sc_num in args.spacecraft_csv.split(","):
    cmd = (
        "python \"" + cmd_path + "\" "
        "--verbose "
        "--url=http://cedar.openmadrigal.org "
        f"--outputDir=./data/{args.run_name}/Satellite_F{sc_num} --user_fullname=ScienceUser "
        "--user_email=noemail@gmail.com "
        '--user_affiliation=NASA '
        '--format=hdf5 '
        f'--startDate=' + args.start_time + ' '
        f'--endDate=' + args.end_time + ' '
        "--inst=8100 "
        f"--kindat=102{sc_num}"
    )

    cprint(f"Downloading data for spacecraft {sc_num} with the following command: ", "green")
    cprint(f"  {cmd}", "cyan")
    subprocess.check_call(cmd, shell=True)

