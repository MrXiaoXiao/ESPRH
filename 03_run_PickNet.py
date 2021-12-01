import numpy as np
from pathlib import Path
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import obspy
import bisect
from pyproj import Geod
import sys
sys.path.append('./src/S_EqT_codes/src/')
sys.path.append('./src/PickNet')

import tensorflow as tf
from data_preprocessing import build_phase_dict
import keras.backend as K
import keras
keras.backend.set_floatx('float32')
import yaml
from random import shuffle
import os
import argparse
from fcn.tester import FCNTester
from scipy import signal

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='03_run_PickNet')
    parser.add_argument('--config-file', dest='config_file', 
                        type=str, help='Configuration file path',default='./default_pipline_config.yaml')
    args = parser.parse_args()
    cfgs = yaml.load(open(args.config_file,'r'),Loader=yaml.SafeLoader)
    task_dir = './' + cfgs['TASKID'] + '/'
    os.chdir(task_dir)
    import tensorflow as tf
    import keras.backend.tensorflow_backend as KTF
    def get_session(gpu_fraction=0.4):
        gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=gpu_fraction)
        return tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))
    session = get_session()
    KTF.set_session(session)
    #os.environ['CUDA_VISIBLE_DEVICES']  = cfgs['EqT']['gpuid']

    #session = tf.Session(config=tf.ConfigProto())
    tester = FCNTester('../src/PickNet/fcn/configs/test_config_example_P_pipeline.yaml')
    tester.setup(session)
    tester.setup_testing_for_pipeline(session)

    wavetype = 'P'
    print('On P-- {}'.format(cfgs['TASKID']))
    reftime = obspy.UTCDateTime(cfgs['REAL']['ref_time'])
    sta_file_name = cfgs['REAL']['save_sta']
    raw_folder_name = './' + cfgs['EqT']['mseed_dir']
    res_folder_name = cfgs['S_EqT']['txt_folder']
    
    if os.path.exists(cfgs['PickNet']['txt_folder']):
        pass
    else:
        os.makedirs(cfgs['PickNet']['txt_folder'])

    phase_dict, station_list = build_phase_dict(sta_file_name, res_folder_name, wavetype='P')
    for sta_info in station_list:
        print('On {}'.format(sta_info))
        sta_name = sta_info[1][3:]
        sta_path = Path(raw_folder_name + '/{}/'.format(sta_name))
        t_waveform = None
        t_wavekey = None

        f_picknet = open('./{}/{}.{}.txt'.format(cfgs['PickNet']['txt_folder'],sta_info[1],wavetype),'w')

        for phase_dx, phase_time in enumerate(phase_dict[sta_info[1]]['{}'.format(wavetype)]):
            seek_time = reftime + phase_time
            #print(seek_time)
            y_key = '{}'.format(seek_time.year)
            m_key = '{}'.format(seek_time.month).rjust(2,'0')
            d_key = '{}'.format(seek_time.day).rjust(2,'0')
            seek_key = '*{}*T*T*.mseed'.format(y_key+m_key+d_key)
            if t_wavekey == seek_key:
                pass
            else:
                t_wavekey = seek_key
                st = obspy.Stream()
                for waveform in sta_path.glob(seek_key):
                    st += obspy.read(str(waveform))
                st.merge()
                # st.filter('bandpass',freqmin=0.01,freqmax=20,zerophase=True)
                st.sort(keys=['channel'])
            #print(st)
            # sort order ENZ 012
            if len(st.traces) == 3:
                t_st_z = st[2].slice(seek_time - 6, seek_time + 6)
                t_st_z.resample(100.0)
                # perform PickNet - P wave
                if len(t_st_z.data[:]) < 1200:
                    refine_time = seek_time
                else:
                    input_mat = np.zeros([1,1,1200,1])
                    input_mat[0,0,:,0] = t_st_z.data[:1200]
                    input_mat -= np.mean(input_mat)
                    input_mat /= np.max(input_mat)
                    res = tester.run_on_one(input_mat, session)
                    shift_t = np.argmax(res[-1][0,0,540:660,0]) + 540
                    
                    if shift_t <= 550 or shift_t >= 650:
                        refine_time = seek_time
                    else:
                        refine_time = seek_time - 6 + 0.01*shift_t
            else:
                refine_time = seek_time
                # calculate max amp
            write_time = refine_time - reftime
            write_line = '{:.3f} {:.5f} {:.8f}\n'.format(write_time, 1.0, 0.0)

            f_picknet.write(write_line)
        f_picknet.close()
    
    K.clear_session()
    session = get_session()
    KTF.set_session(session)
    #session = tf.Session(config=tf.ConfigProto())
    tester = FCNTester('../src/PickNet/fcn/configs/test_config_example_S_pipeline.yaml')
    tester.setup(session)
    tester.setup_testing_for_pipeline(session)
    wavetype = 'S'
   

    print('On S-- {}'.format(cfgs['TASKID']))
    reftime = obspy.UTCDateTime(cfgs['REAL']['ref_time'])
    sta_file_name = cfgs['REAL']['save_sta']
    raw_folder_name = './' + cfgs['EqT']['mseed_dir']
    res_folder_name = cfgs['S_EqT']['txt_folder']

    if os.path.exists(cfgs['PickNet']['txt_folder']):
        pass
    else:
        os.makedirs(cfgs['PickNet']['txt_folder'])

    phase_dict, station_list = build_phase_dict(sta_file_name, res_folder_name, wavetype)
    for sta_info in station_list:
        print('On {}'.format(sta_info))
        sta_name = sta_info[1][3:]
        sta_path = Path(raw_folder_name + '/{}/'.format(sta_name))
        t_waveform = None
        t_wavekey = None

        f_picknet = open('./{}/{}.{}.txt'.format(cfgs['PickNet']['txt_folder'],sta_info[1],wavetype),'w')

        for phase_dx, phase_time in enumerate(phase_dict[sta_info[1]]['{}'.format(wavetype)]):
            seek_time = reftime + phase_time
            #print(seek_time)
            y_key = '{}'.format(seek_time.year)
            m_key = '{}'.format(seek_time.month).rjust(2,'0')
            d_key = '{}'.format(seek_time.day).rjust(2,'0')
            seek_key = '*{}*T*T*.mseed'.format(y_key+m_key+d_key)
            if t_wavekey == seek_key:
                pass
            else:
                t_wavekey = seek_key
                st = obspy.Stream()
                for waveform in sta_path.glob(seek_key):
                    st += obspy.read(str(waveform))
                st.merge()
                # st.filter('bandpass',freqmin=0.01,freqmax=20,zerophase=True)
                st.sort(keys=['channel'])
            #print(st)
            # sort order ENZ 012
            if len(st.traces) == 3:
                t_st_e = st[0].slice(seek_time - 8, seek_time + 8)
                t_st_e.resample(100.0)
                t_st_n = st[1].slice(seek_time - 8, seek_time + 8)
                t_st_n.resample(100.0)                    
                # perform PickNet - S wave
                # print(len(t_st_e.data[:]))
                if len(t_st_e.data[:]) < 1600 or len(t_st_n.data[:]) < 1600: 
                    refine_time = seek_time
                else:
                    input_mat = np.zeros([1,1,1600,2])

                    input_mat[0,0,:,0] = t_st_e.data[:1600]
                    input_mat[0,0,:,0] -= np.mean(input_mat[0,0,:,0])
                    input_mat[0,0,:,0] /= np.max(input_mat[0,0,:,0])

                    input_mat[0,0,:,1] = t_st_n.data[:1600]
                    input_mat[0,0,:,1] -= np.mean(input_mat[0,0,:,1])
                    input_mat[0,0,:,1] /= np.max(input_mat[0,0,:,1])

                    res = tester.run_on_one(input_mat, session)
                    shift_t = np.argmax(res[-1][0,0,690:911,0]) + 690
                    #print(shift_t)
                    if shift_t <= 700 or shift_t >= 900:
                        refine_time = seek_time
                    else:
                        refine_time = seek_time - 8 + 0.01*shift_t
            else:
                refine_time = seek_time
            write_time = refine_time - reftime
            write_line = '{:.3f} {:.5f} {:.8f}\n'.format(write_time, 1.0, 0.0)
            f_picknet.write(write_line)
        f_picknet.close()