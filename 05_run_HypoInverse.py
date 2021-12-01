import obspy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
from obspy import read, UTCDateTime
import os
import math
import yaml
from pathlib import Path
import shutil
import re
import argparse

def create_station_file(cfgs):
    real_sta_file_path = cfgs['REAL']['save_sta']
    real_sta_file = open(real_sta_file_path, 'r')
    
    hypo_sta_file_path = cfgs['HypoInverse']['save_sta']
    hypo_sta_file = open(hypo_sta_file_path, 'w')

    for line in real_sta_file.readlines():
        splits = line.split(' ')
        lat = '{:.5f}'.format(float(splits[1]))
        lon = '{:.5f}'.format(float(splits[0]))
        code = splits[2]+'.'+splits[3]
        ele = '{:.2f}'.format(float(splits[-1])*1000.0)
        pad = '-1'

        hypo_line = '\t'.join([code,lat,lon,ele,pad]) + '\n'

        hypo_sta_file.write(hypo_line)
    
    real_sta_file.close()
    hypo_sta_file.close()
    return

def create_pha_file(cfgs):
    real_event_dict_path =  cfgs['REAL']['eqt_catalog_dir'] + cfgs['HypoInverse']['eqt_event_dict']
    real_event_dict = np.load(real_event_dict_path,allow_pickle=True )[()]
    
    hypo_pha_file_path = cfgs['HypoInverse']['save_pha_eqt']
    hypo_pha_file = open(hypo_pha_file_path,'w')
    
    for e_key in real_event_dict.keys():
        ot = str(real_event_dict[e_key]['REAL_TIME'])
        lat = '{:5f}'.format(real_event_dict[e_key]['REAL_LAT'])
        lon = '{:5f}'.format(real_event_dict[e_key]['REAL_LON'])
        dep = '{:5f}'.format(real_event_dict[e_key]['REAL_DEP'])
        mag = '1.0'
        # create event line
        event_line = ','.join([ot,lat,lon,dep,mag]) + '\n'
        hypo_pha_file.write(event_line)
        
        temp_pick_dict = dict()
        for pick_info in real_event_dict[e_key]['Picks']:
            code = pick_info[0]
            pick_type = pick_info[1]
            pick_time = pick_info[2]

            if code in temp_pick_dict.keys():
                temp_pick_dict[code][pick_type] = pick_time
            else:
                temp_pick_dict[code] = dict()
                temp_pick_dict[code]['P'] = -1
                temp_pick_dict[code]['S'] = -1
                temp_pick_dict[code][pick_type] = pick_time
        
        for pick_key in temp_pick_dict.keys():
            net = pick_key.split('.')[0]
            sta = pick_key.split('.')[1]
            tp = str(temp_pick_dict[pick_key]['P'])
            ts = str(temp_pick_dict[pick_key]['S'])
            pick_line =  ','.join([net,sta,tp,ts]) + ',-1,-1,-1\n'
            hypo_pha_file.write(pick_line)
    
    hypo_pha_file.close()

    real_event_dict_path =  cfgs['REAL']['seqt_catalog_dir'] + cfgs['HypoInverse']['seqt_event_dict']
    real_event_dict = np.load(real_event_dict_path,allow_pickle=True )[()]
    
    hypo_pha_file_path = cfgs['HypoInverse']['save_pha_seqt']
    hypo_pha_file = open(hypo_pha_file_path,'w')
    
    for e_key in real_event_dict.keys():
        ot = str(real_event_dict[e_key]['REAL_TIME'])
        lat = '{:5f}'.format(real_event_dict[e_key]['REAL_LAT'])
        lon = '{:5f}'.format(real_event_dict[e_key]['REAL_LON'])
        dep = '{:5f}'.format(real_event_dict[e_key]['REAL_DEP'])
        mag = '1.0'
        # create event line
        event_line = ','.join([ot,lat,lon,dep,mag]) + '\n'
        hypo_pha_file.write(event_line)
        
        temp_pick_dict = dict()
        for pick_info in real_event_dict[e_key]['Picks']:
            code = pick_info[0]
            pick_type = pick_info[1]
            pick_time = pick_info[2]

            if code in temp_pick_dict.keys():
                temp_pick_dict[code][pick_type] = pick_time
            else:
                temp_pick_dict[code] = dict()
                temp_pick_dict[code]['P'] = -1
                temp_pick_dict[code]['S'] = -1
                temp_pick_dict[code][pick_type] = pick_time
        
        for pick_key in temp_pick_dict.keys():
            net = pick_key.split('.')[0]
            sta = pick_key.split('.')[1]
            tp = str(temp_pick_dict[pick_key]['P'])
            ts = str(temp_pick_dict[pick_key]['S'])
            pick_line =  ','.join([net,sta,tp,ts]) + ',-1,-1,-1\n'
            hypo_pha_file.write(pick_line)
    
    hypo_pha_file.close()

    real_event_dict_path =  cfgs['REAL']['picknet_catalog_dir'] + cfgs['HypoInverse']['picknet_event_dict']
    real_event_dict = np.load(real_event_dict_path,allow_pickle=True )[()]
    
    hypo_pha_file_path = cfgs['HypoInverse']['save_pha_picknet']
    hypo_pha_file = open(hypo_pha_file_path,'w')
    
    for e_key in real_event_dict.keys():
        ot = str(real_event_dict[e_key]['REAL_TIME'])
        lat = '{:5f}'.format(real_event_dict[e_key]['REAL_LAT'])
        lon = '{:5f}'.format(real_event_dict[e_key]['REAL_LON'])
        dep = '{:5f}'.format(real_event_dict[e_key]['REAL_DEP'])
        mag = '1.0'
        # create event line
        event_line = ','.join([ot,lat,lon,dep,mag]) + '\n'
        hypo_pha_file.write(event_line)
        
        temp_pick_dict = dict()
        for pick_info in real_event_dict[e_key]['Picks']:
            code = pick_info[0]
            pick_type = pick_info[1]
            pick_time = pick_info[2]

            if code in temp_pick_dict.keys():
                temp_pick_dict[code][pick_type] = pick_time
            else:
                temp_pick_dict[code] = dict()
                temp_pick_dict[code]['P'] = -1
                temp_pick_dict[code]['S'] = -1
                temp_pick_dict[code][pick_type] = pick_time
        
        for pick_key in temp_pick_dict.keys():
            net = pick_key.split('.')[0]
            sta = pick_key.split('.')[1]
            tp = str(temp_pick_dict[pick_key]['P'])
            ts = str(temp_pick_dict[pick_key]['S'])
            pick_line =  ','.join([net,sta,tp,ts]) + ',-1,-1,-1\n'
            hypo_pha_file.write(pick_line)
    
    hypo_pha_file.close()

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='00_data_download')
    parser.add_argument('--config-file', dest='config_file', 
                        type=str, help='Configuration file path',default='./default_pipline_config.yaml')
    args = parser.parse_args()
    cfgs = yaml.load(open(args.config_file,'r'),Loader=yaml.SafeLoader)

    task_dir = './' + cfgs['TASKID'] + '/'
    os.chdir(task_dir)

    create_station_file(cfgs)
    create_pha_file(cfgs)

    os.rename(cfgs['HypoInverse']['save_sta'], '../HypoInverse_scripts/input/HYPO.sta')
    os.rename(cfgs['HypoInverse']['save_pha_eqt'], '../HypoInverse_scripts/input/HYPO.pha')
    os.chdir('../HypoInverse_scripts')
    hypo_output = os.system('python run_hyp.py')
    print('STATUS: {}'.format(hypo_output))
    os.chdir('..')
    os.chdir(task_dir)
    os.rename( '../HypoInverse_scripts/output/example.ctlg', '{}/eqt_hypoInverse.ctlg'.format(cfgs['REAL']['eqt_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example.pha','{}/eqt_hypoInverse.pha'.format(cfgs['REAL']['eqt_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example.sum','{}/eqt_hypoInverse.sum'.format(cfgs['REAL']['eqt_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example_good.csv','{}/eqt_hypoInverse.good'.format(cfgs['REAL']['eqt_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example_bad.csv','./{}/eqt_hypoInverse.bad'.format(cfgs['REAL']['eqt_catalog_dir']))

    os.rename(cfgs['HypoInverse']['save_pha_seqt'], '../HypoInverse_scripts/input/HYPO.pha')
    os.chdir('../HypoInverse_scripts')
    hypo_output = os.system('python run_hyp.py')
    print('STATUS: {}'.format(hypo_output))
    os.chdir('..')
    os.chdir(task_dir)
    os.rename( '../HypoInverse_scripts/output/example.ctlg', '{}/seqt_hypoInverse.ctlg'.format(cfgs['REAL']['seqt_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example.pha','{}/seqt_hypoInverse.pha'.format(cfgs['REAL']['seqt_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example.sum','{}/seqt_hypoInverse.sum'.format(cfgs['REAL']['seqt_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example_good.csv','{}/seqt_hypoInverse.good'.format(cfgs['REAL']['seqt_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example_bad.csv','{}/seqt_hypoInverse.bad'.format(cfgs['REAL']['seqt_catalog_dir']))

    os.rename(cfgs['HypoInverse']['save_pha_picknet'], '../HypoInverse_scripts/input/HYPO.pha')
    os.chdir('../HypoInverse_scripts')
    hypo_output = os.system('python run_hyp.py')
    print('STATUS: {}'.format(hypo_output))
    os.chdir('..')
    os.chdir(task_dir)
    os.rename( '../HypoInverse_scripts/output/example.ctlg', '{}/picknet_hypoInverse.ctlg'.format(cfgs['REAL']['picknet_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example.pha','{}/picknet_hypoInverse.pha'.format(cfgs['REAL']['picknet_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example.sum','{}/picknet_hypoInverse.sum'.format(cfgs['REAL']['picknet_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example_good.csv','{}/picknet_hypoInverse.good'.format(cfgs['REAL']['picknet_catalog_dir']))
    os.rename( '../HypoInverse_scripts/output/example_bad.csv','{}/picknet_hypoInverse.bad'.format(cfgs['REAL']['picknet_catalog_dir']))