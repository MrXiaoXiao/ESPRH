import os
import sys
import argparse
import yaml

import time
import numpy as np

import tensorflow as tf
import math
from fcn.models.vgg16_hed import Vgg16_HED
from fcn.models.vgg16_srn import Vgg16_SRN
from fcn.models.picknet import PickNet
from fcn.utils.io import IO

import matplotlib.pyplot as plt
plt.switch_backend('agg')

class FCNTester():

    def __init__(self, config_file):

        self.io = IO()
        self.init = True

        try:
            pfile = open(config_file)
            self.cfgs = yaml.load(pfile, Loader=yaml.SafeLoader)
            pfile.close()

        except Exception as err:

            self.io.print_error('Error reading config file {}, {}'.format(config_file), err)

    def setup(self, session):

        try:
            if self.cfgs['using_model'] == 'HED':               
                self.model = Vgg16_HED(self.cfgs, run='testing')
                self.io.print_info('Done initializing VGG-16_HED model')
            elif self.cfgs['using_model'] == 'SRN':               
                self.model = Vgg16_SRN(self.cfgs, run='testing')
                self.io.print_info('Done initializing VGG-16_SRN model')
            elif self.cfgs['using_model'] == 'PickNet':
                self.model = PickNet(self.cfgs, run='testing')
                self.io.print_info('Done initializing PickNet model')
            else:
                self.io.print_info('Fail to setup model')
                return

            meta_model_file = os.path.join(self.cfgs['save_dir'], 'models/{}-model-{}'.format(self.cfgs['using_model'],self.cfgs['test_snapshot']))

            saver = tf.train.Saver()
            saver.restore(session, meta_model_file)

            self.io.print_info('Done restoring model from {}'.format(meta_model_file))

        except Exception as err:
            print('WTF?')
            self.io.print_error('Error setting up model, {}'.format(err))
            self.init = False
    
    def setup_testing_for_pipeline(self, session):
        if not self.init:
            print('No INIT')
            return
        self.model.setup_testing(session)

        return
    
    def run_on_one(self, input_data, session):
        outputs = session.run(self.model.predictions, feed_dict={self.model.inputs: input_data})
        return outputs
    
    def run(self, session):

        if not self.init:
            print('No INIT')
            return
        
        self.model.setup_testing(session)

        test_file = self.cfgs['testing']['filename'] + '.npy'
        test_threshold = float(self.cfgs['testing_threshold'])
        self.testing_data = np.load(test_file)
        
        self.n_samples = len(self.testing_data)
        
        self.testing_data = np.reshape(self.testing_data,(self.n_samples,self.cfgs['testing']['data_width'],self.cfgs['testing']['n_channels']))
        
        print(self.testing_data.shape)
        
        batch_size = int(self.cfgs['batch_size_test'])
        total_step = int(math.ceil(self.n_samples/batch_size))        
        
        self.io.print_info('Start Testing')
        
        total_output = list()
        
        fuse_picks = list()
        t0 = time.clock()
        for i in range(total_step):
            self.io.print_info('Step: {} of {}'.format(i,total_step))
            
            i_start = i * batch_size
            i_end = (i+1) * batch_size
            
            if i_end > self.n_samples:
                im = self.testing_data[i_start: self.n_samples,:,:]
                im = np.pad(im,pad_width=((0,i_end-self.n_samples),(0,0),(0,0)),mode='constant')
                #i_end = self.n_samples
            else:
                im = self.testing_data[i_start:i_end,:,:]
            
            
            im = np.reshape(im,
                            (i_end-i_start,
                             self.cfgs['testing']['data_height'],
                             self.cfgs['testing']['data_width'],
                             self.cfgs['testing']['n_channels']))
                       
            outputs = session.run(self.model.predictions, feed_dict={self.model.inputs: im})
            
            total_output.append(outputs[-1])
            #print(im.shape)
            #print(outputs[5].shape)
            
            for j in range(i_end-i_start):
                if outputs[-1][j,0,self.cfgs['testing']['limit_left']:self.cfgs['testing']['limit_right'],0].max() < test_threshold:
                    continue
                pick_point_fuse = self.cfgs['testing']['limit_left'] + outputs[-1][j,0,self.cfgs['testing']['limit_left']:self.cfgs['testing']['limit_right'],0].argmax()                       
                
                fuse_picks.append('{}_{}'.format(i*batch_size + j,pick_point_fuse))
        self.io.print_info('{} seconds'.format(time.clock() - t0))
        np.save(self.cfgs['testing']['filename']+'_output.npy',total_output)
        np.save(self.cfgs['testing']['filename']+'_fuse_picks.npy',fuse_picks)
        self.io.print_info('Done testing {}'.format(test_file))
        self.io.print_info('Pick {} out of {}'.format(len(fuse_picks),self.n_samples))
        
