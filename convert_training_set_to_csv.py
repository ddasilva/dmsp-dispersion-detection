#!/usr/bin/python

from datetime import date, timedelta
import os
import pandas as pd
from sklearn.model_selection import train_test_split
import subprocess
import progressbar

with_dir = './Examples/WithDispersion'
without_dir = './Examples/WithoutDispersion'

with_files = [f for f in os.listdir(with_dir) if f.endswith('.png')]
without_files = [f for f in os.listdir(without_dir) if f.endswith('.png')]

rows = []


for cls, file_list in [(0, without_files), (1, with_files)]:
    bar = progressbar.ProgressBar()

    for file in bar(file_list):
        toks = file.split('_')        
        
        toks2_toks = toks[2].replace('.png', '').split('-')
        sat = int(toks[0][1:3])
        start_hour = int(toks2_toks[0][:2])
        start_min = int(toks2_toks[0][2:])
        end_hour = int(toks2_toks[1][:2])
        end_min = int(toks2_toks[1][2:])
        start_year = int(toks[1][:4])
        start_mon = int(toks[1][4:6])
        start_day = int(toks[1][6:8])        
        
        prev_day = date(start_year, start_mon, start_day) - timedelta(days=1)


        if start_hour >= 60:
            continue
        
        filename = f'data/{start_year:04d}/{start_mon:02d}/{start_day:02d}/dms_{start_year:04d}{start_mon:02d}{start_day:02d}_{sat:02d}e.001.hdf5'
        prev_filename = f'data/{start_year:04d}/{start_mon:02d}/{start_day:02d}/dms_{prev_day.year:04d}{prev_day.month:02d}{prev_day.day:02d}_{sat:02d}e.001.hdf5'

        if not os.path.exists(filename):
            # --kindat documentation available at http://cedar.openmadrigal.org/kindatMetadata/
            cmd = f"/home/daniel/anaconda3/bin/globalDownload.py --verbose --url=http://cedar.openmadrigal.org --outputDir=./data/{start_year:04d}/{start_mon:02d}/{start_day:02d} --user_fullname='Daniel+da+Silva' --user_email=daniel.e.dasilva@nasa.gov --user_affiliation='NASA' --format='hdf5' --startDate='{start_mon:02d}/{start_day:02d}/{start_year:04d}' --endDate='{start_mon:02d}/{start_day:02d}/{start_year:04d}' --inst=8100 --kindat=102{sat}"

            subprocess.check_call(cmd, shell=True)

        if os.path.exists(prev_filename):
            print("Removing previous day file")
            os.remove(prev_filename)

        rows.append((
            cls,
            filename,
            int(sat),
            f'{start_year:04d}-{start_mon:02d}-{start_day:02d}',
            f'{start_year:04d}-{start_mon:02d}-{start_day:02d}T{start_hour:02d}:{start_min:02d}Z',
            f'{start_year:04d}-{start_mon:02d}-{start_day:02d}T{end_hour:02d}:{end_min:02d}Z'
        ))

                
df = pd.DataFrame(rows,
                  columns=['class', 'filename', 'sat', 'date', 'start_time', 'end_time'])

train, val = train_test_split(df, test_size=0.1, random_state=42, shuffle=True)

train.to_csv('data/train.csv', index=False)
val.to_csv('data/val.csv', index=False)




    
