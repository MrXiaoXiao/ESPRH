# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 19:17:05 2018
"""
import os
import sys
import yaml
import argparse
import tensorflow as tf
import numpy as np

from fcn.utils.io import IO
from fcn.data.data_parser import DataParser
from fcn.models.vgg16_hed import Vgg16_HED
from fcn.models.vgg16_srn import Vgg16_SRN
from fcn.models.picknet import PickNet

class FCNTrainer():
    def __init__(self,config_file):
        #set IO
        self.io = IO()
        #set init status
        self.init = True
        
        #open config file
        try:
            pfile = open(config_file)
            self.cfgs = yaml.load(pfile)
            pfile.close()
        
        except Exception as err:
            print('Error reading config file {}, {}'.format(config_file, err))
            
    def setup(self):
        #set up inference model
        try:
            if self.cfgs['using_model'] == 'HED':                
                #init vgg16 model
                self.model = Vgg16_HED(self.cfgs)
                self.io.print_info('Done initializing VGG-16_HED model')
            elif self.cfgs['using_model'] == 'SRN':
                self.model = Vgg16_SRN(self.cfgs)
                self.io.print_info('Done initializing VGG-16_SRN model')
            elif self.cfgs['using_model'] == 'PickNet':
                self.model = PickNet(self.cfgs)
                self.io.print_info('Done initializing PickNet model')
            else:
                self.io.print_info('Fail to setup model')
                return
        
        except Exception as err:           
            print('Error setting up inference model, {}'.format(err))
            
    def run(self, session):
        loss_list = list()
        
        if not self.init:
            return
        #load train data
        train_data = DataParser(self.cfgs)
        
        #start training
        self.model.setup_training(session)

        opt = tf.train.AdamOptimizer(self.cfgs['optimizer_params']['learning_rate'])
        train = opt.minimize(self.model.loss)

        session.run(tf.global_variables_initializer())

        for idx in range(self.cfgs['max_iterations']):
            
            waves, picks = train_data.get_training_batch()

            _, summary, loss = session.run([train, self.model.merged_summary, self.model.loss],
                                           feed_dict={self.model.inputs: waves, self.model.groundtruths: picks})

            self.model.train_writer.add_summary(summary, idx)
            
            
            if idx % self.cfgs['print_interval'] == 0:    
                self.io.print_info('[{}/{}] TRAINING loss : {}'.format(idx, self.cfgs['max_iterations'], loss))

            if idx % self.cfgs['save_interval'] == 0:
                saver = tf.train.Saver()
                saver.save(session, os.path.join(self.cfgs['save_dir'], 'models/{}-model'.format(self.cfgs['using_model'])), global_step=idx)
                loss_list.append(loss)
                               
            if idx % self.cfgs['val_interval'] == 0:
                waves, picks = train_data.get_validation_batch()
                summary, error = session.run([self.model.merged_summary, self.model.error], feed_dict={self.model.inputs: waves, self.model.groundtruths: picks})
                self.model.val_writer.add_summary(summary, idx)
                self.io.print_info('[{}/{}] VALIDATION error : {}'.format(idx, self.cfgs['max_iterations'], error))
            if idx % 5000 == 0:
                np.save(self.cfgs['save_dir']+'/losses.npy',loss_list)
        self.model.train_writer.close()
        np.save(self.cfgs['save_dir']+'/losses.npy',loss_list)
