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

"""
Workaround codes. Need Cleaning.
"""

def convert2sec(t, t_ref):
    """
    convert UTCDatetime object to seconds
    Params:
    t       UTCDateTime     Time to be converted
    t_ref   UTCDateTime     Reference time
    """
    t_utc = UTCDateTime(t)
    t_ref_utc = UTCDateTime(t_ref)
    return t_utc - t_ref_utc

def runREAL(cfgs):
    """
    Run REAL Scripts
    """
    if os.path.exists(cfgs['REAL']['eqt_catalog_dir']):
        pass
    else:
        os.makedirs(cfgs['REAL']['eqt_catalog_dir'])

    if os.path.exists(cfgs['REAL']['seqt_catalog_dir']):
        pass
    else:
        os.makedirs(cfgs['REAL']['seqt_catalog_dir'])

    if os.path.exists(cfgs['REAL']['picknet_catalog_dir']):
        pass
    else:
        os.makedirs(cfgs['REAL']['picknet_catalog_dir'])

    for idx in range(len(cfgs['REAL']['year'])):
        # copy temp perl file
        f_perl = open('../REAL_scripts/runREAL.pl', 'r')
        f_perl_source = f_perl.read()
        f_perl.close()
        f_perl_source = f_perl_source.replace('YEAR_KEY', cfgs['REAL']['year'][idx])
        f_perl_source = f_perl_source.replace('MON_KEY', cfgs['REAL']['mon'][idx])
        f_perl_source = f_perl_source.replace('DAY_KEY', cfgs['REAL']['day'][idx])
        f_perl_source = f_perl_source.replace('DIR_KEY','\"'  + cfgs['REAL']['eqt_dir']  + '\"')
        f_perl_source = f_perl_source.replace('STATION_KEY', cfgs['REAL']['station'])
        f_perl_source = f_perl_source.replace('TTIME_KEY', cfgs['REAL']['ttime'])
        f_perl_source = f_perl_source.replace('R_KEY', cfgs['REAL']['R'])
        f_perl_source = f_perl_source.replace('G_KEY', cfgs['REAL']['G'])
        f_perl_source = f_perl_source.replace('V_KEY', cfgs['REAL']['V'])
        f_perl_source = f_perl_source.replace('S_KEY', cfgs['REAL']['S'])
        f_perl_temp = open('../REAL_scripts/runREAL_temp.pl','w')
        f_perl_temp.write(f_perl_source)
        f_perl_temp.close()
        real_output = os.system('../REAL_scripts/runREAL_temp.pl')
        print('STATUS: {}'.format(real_output))
        
        os.rename('./catalog_sel.txt', '{}eqt_real_catalog_sel.txt'.format(cfgs['REAL']['eqt_catalog_dir']))
        os.rename('./phase_sel.txt', '{}eqt_real_phase_sel.txt'.format(cfgs['REAL']['eqt_catalog_dir']))

        # copy temp perl file
        f_perl = open('../REAL_scripts/runREAL.pl', 'r')
        f_perl_source = f_perl.read()
        f_perl.close()
        f_perl_source = f_perl_source.replace('YEAR_KEY', cfgs['REAL']['year'][idx])
        f_perl_source = f_perl_source.replace('MON_KEY', cfgs['REAL']['mon'][idx])
        f_perl_source = f_perl_source.replace('DAY_KEY', cfgs['REAL']['day'][idx])
        f_perl_source = f_perl_source.replace('DIR_KEY','\"'  + cfgs['REAL']['seqt_dir']  + '\"')
        f_perl_source = f_perl_source.replace('STATION_KEY', cfgs['REAL']['station'])
        f_perl_source = f_perl_source.replace('TTIME_KEY', cfgs['REAL']['ttime'])
        f_perl_source = f_perl_source.replace('R_KEY', cfgs['REAL']['R'])
        f_perl_source = f_perl_source.replace('G_KEY', cfgs['REAL']['G'])
        f_perl_source = f_perl_source.replace('V_KEY', cfgs['REAL']['V'])
        f_perl_source = f_perl_source.replace('S_KEY', cfgs['REAL']['S'])
        f_perl_temp = open('../REAL_scripts/runREAL_temp.pl','w')
        f_perl_temp.write(f_perl_source)
        f_perl_temp.close()
        real_output = os.system('../REAL_scripts/runREAL_temp.pl')
        print('STATUS: {}'.format(real_output))
        
        os.rename('./catalog_sel.txt', '{}seqt_real_catalog_sel.txt'.format(cfgs['REAL']['seqt_catalog_dir']))
        os.rename('./phase_sel.txt', '{}seqt_real_phase_sel.txt'.format(cfgs['REAL']['seqt_catalog_dir']))

        # copy temp perl file
        f_perl = open('../REAL_scripts/runREAL.pl', 'r')
        f_perl_source = f_perl.read()
        f_perl.close()
        f_perl_source = f_perl_source.replace('YEAR_KEY', cfgs['REAL']['year'][idx])
        f_perl_source = f_perl_source.replace('MON_KEY', cfgs['REAL']['mon'][idx])
        f_perl_source = f_perl_source.replace('DAY_KEY', cfgs['REAL']['day'][idx])
        f_perl_source = f_perl_source.replace('DIR_KEY','\"'  + cfgs['REAL']['picknet_dir']  + '\"')
        f_perl_source = f_perl_source.replace('STATION_KEY', cfgs['REAL']['station'])
        f_perl_source = f_perl_source.replace('TTIME_KEY', cfgs['REAL']['ttime'])
        f_perl_source = f_perl_source.replace('R_KEY', cfgs['REAL']['R'])
        f_perl_source = f_perl_source.replace('G_KEY', cfgs['REAL']['G'])
        f_perl_source = f_perl_source.replace('V_KEY', cfgs['REAL']['V'])
        f_perl_source = f_perl_source.replace('S_KEY', cfgs['REAL']['S'])
        f_perl_temp = open('../REAL_scripts/runREAL_temp.pl','w')
        f_perl_temp.write(f_perl_source)
        f_perl_temp.close()
        real_output = os.system('../REAL_scripts/runREAL_temp.pl')
        print('STATUS: {}'.format(real_output))
        
        os.rename('./catalog_sel.txt', '{}picknet_real_catalog_sel.txt'.format(cfgs['REAL']['picknet_catalog_dir']))
        os.rename('./phase_sel.txt', '{}picknet_real_phase_sel.txt'.format(cfgs['REAL']['picknet_catalog_dir']))

    return

def merge_phasesel(cfgs):
    """
    Merge phase sel files
    """
    e_dict = dict()
    base_time = obspy.UTCDateTime(cfgs['REAL']['ref_time'])
    e_ID = None
    f_sel = open('{}/eqt_real_phase_sel.txt'.format(cfgs['REAL']['eqt_catalog_dir']),'r')
    for line in f_sel.readlines():
        line_split = re.sub('\s{2,}',' ',line).split(' ')
        if len(line_split) > 11:
            e_ID = '{}'.format(int(line_split[1]))
            e_dict[e_ID] = dict()
            real_time = base_time + float(line_split[6])
            e_dict[e_ID]['REAL_TIME'] = real_time
            e_dict[e_ID]['REAL_LAT'] = float(line_split[8])
            e_dict[e_ID]['REAL_LON'] = float(line_split[9])
            e_dict[e_ID]['REAL_DEP'] = float(line_split[10])
            e_dict[e_ID]['Picks'] = list()
        else:
            sta_name = line_split[1] + '.' + line_split[2]
            pick_type = line_split[3]
            pick_time = base_time + float(line_split[4])
            if pick_time - e_dict[e_ID]['REAL_TIME'] < 0.01:
                continue
            e_dict[e_ID]['Picks'].append([sta_name, pick_type, pick_time])
    f_sel.close()
    #print_dict(e_dict)
    np.save('{}/eqt_real_e_dict.npy'.format(cfgs['REAL']['eqt_catalog_dir']),e_dict)

    e_dict = dict()
    base_time = obspy.UTCDateTime(cfgs['REAL']['ref_time'])
    e_ID = None
    f_sel = open('{}/seqt_real_phase_sel.txt'.format(cfgs['REAL']['seqt_catalog_dir']),'r')
    for line in f_sel.readlines():
        line_split = re.sub('\s{2,}',' ',line).split(' ')
        if len(line_split) > 11:
            e_ID = '{}'.format(int(line_split[1]))
            e_dict[e_ID] = dict()
            real_time = base_time + float(line_split[6])
            e_dict[e_ID]['REAL_TIME'] = real_time
            e_dict[e_ID]['REAL_LAT'] = float(line_split[8])
            e_dict[e_ID]['REAL_LON'] = float(line_split[9])
            e_dict[e_ID]['REAL_DEP'] = float(line_split[10])
            e_dict[e_ID]['Picks'] = list()
        else:
            sta_name = line_split[1] + '.' + line_split[2]
            pick_type = line_split[3]
            pick_time = base_time + float(line_split[4])
            if pick_time - e_dict[e_ID]['REAL_TIME'] < 0.01:
                continue
            e_dict[e_ID]['Picks'].append([sta_name, pick_type, pick_time])
    f_sel.close()
    
    np.save('{}/seqt_real_e_dict.npy'.format(cfgs['REAL']['seqt_catalog_dir']),e_dict)

    e_dict = dict()
    base_time = obspy.UTCDateTime(cfgs['REAL']['ref_time'])
    e_ID = None
    f_sel = open('{}/picknet_real_phase_sel.txt'.format(cfgs['REAL']['picknet_catalog_dir']),'r')
    for line in f_sel.readlines():
        line_split = re.sub('\s{2,}',' ',line).split(' ')
        if len(line_split) > 11:
            e_ID = '{}'.format(int(line_split[1]))
            e_dict[e_ID] = dict()
            real_time = base_time + float(line_split[6])
            e_dict[e_ID]['REAL_TIME'] = real_time
            e_dict[e_ID]['REAL_LAT'] = float(line_split[8])
            e_dict[e_ID]['REAL_LON'] = float(line_split[9])
            e_dict[e_ID]['REAL_DEP'] = float(line_split[10])
            e_dict[e_ID]['Picks'] = list()
        else:
            sta_name = line_split[1] + '.' + line_split[2]
            pick_type = line_split[3]
            pick_time = base_time + float(line_split[4])
            if pick_time - e_dict[e_ID]['REAL_TIME'] < 0.01:
                continue
            e_dict[e_ID]['Picks'].append([sta_name, pick_type, pick_time])
    f_sel.close()
    
    np.save('{}/picknet_real_e_dict.npy'.format(cfgs['REAL']['picknet_catalog_dir']),e_dict)
    return

def print_dict(e_dict):
    for key in e_dict.keys():
        print('E_ID: {} VELTime: {} LAT: {} LON: {} DEP: {}'.format(key,
                                                                e_dict[key]['VELEST_TIME'],
                                                                e_dict[key]['VELEST_LAT'],
                                                                e_dict[key]['VELEST_LON'],
                                                                e_dict[key]['VELEST_DEP']))

        print('REALTime: {} LAT: {} LON: {} DEP: {}'.format(e_dict[key]['REAL_TIME'],
                                                                e_dict[key]['REAL_LAT'],
                                                                e_dict[key]['REAL_LON'],
                                                                e_dict[key]['REAL_DEP']))
        for pick in e_dict[key]['Picks']:
            print(pick)
    return

def pad_empty_sta(cfgs):
    f = open(cfgs['REAL']['save_sta'],'r')
    lines = f.readlines()
    f.close()

    save_folder = cfgs['REAL']['eqt_dir']
    for line in lines:
        splits = line.split(' ')
        sta_name = splits[3]
        net_name = splits[2]

        t_P_name = save_folder + net_name + '.' +sta_name+'.P.txt'
        t_S_name = save_folder + net_name + '.' +sta_name+'.S.txt'
        if os.path.exists(t_P_name):
            pass
        else:
            t_f = open(t_P_name, 'w')
            t_f.close()
        if os.path.exists(t_S_name):
            pass
        else:
            t_f = open(t_S_name, 'w')
            t_f.close()

    f = open(cfgs['REAL']['save_sta'],'r')
    lines = f.readlines()
    f.close()    
    save_folder = cfgs['REAL']['seqt_dir']
    for line in lines:
        splits = line.split(' ')
        sta_name = splits[3]
        net_name = splits[2]

        t_P_name = save_folder + net_name + '.' +sta_name+'.P.txt'
        t_S_name = save_folder + net_name + '.' +sta_name+'.S.txt'
        if os.path.exists(t_P_name):
            pass
        else:
            t_f = open(t_P_name, 'w')
            t_f.close()
        if os.path.exists(t_S_name):
            pass
        else:
            t_f = open(t_S_name, 'w')
            t_f.close()
    
    f = open(cfgs['REAL']['save_sta'],'r')
    lines = f.readlines()
    f.close()    
    save_folder = cfgs['REAL']['picknet_dir']
    for line in lines:
        splits = line.split(' ')
        sta_name = splits[3]
        net_name = splits[2]

        t_P_name = save_folder + net_name + '.' +sta_name+'.P.txt'
        t_S_name = save_folder + net_name + '.' +sta_name+'.S.txt'
        if os.path.exists(t_P_name):
            pass
        else:
            t_f = open(t_P_name, 'w')
            t_f.close()
        if os.path.exists(t_S_name):
            pass
        else:
            t_f = open(t_S_name, 'w')
            t_f.close()
    return
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='04_run_REAL')
    parser.add_argument('--config-file', dest='config_file', 
                        type=str, help='Configuration file path',default='./default_pipline_config.yaml')
    args = parser.parse_args()
    cfgs = yaml.load(open(args.config_file,'r'),Loader=yaml.SafeLoader)
    task_dir = './' + cfgs['TASKID'] + '/'
    os.chdir(task_dir)
    pad_empty_sta(cfgs)
    runREAL(cfgs)
    merge_phasesel(cfgs)
   
    