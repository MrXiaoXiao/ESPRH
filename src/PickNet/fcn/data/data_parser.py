import os
import sys
import time
import numpy as np
from fcn.utils.io import IO

class DataParser():

    def __init__(self, cfgs, run='training'):

        self.io = IO()
        self.cfgs = cfgs
        self.is_test = False
        self.data_file_name = cfgs[run]['filename']+'.npy'

        self.io.print_info('Training or Testing data set-up from {}'.format(os.path.join(self.data_file_name)))
        self.total_data = np.load(self.data_file_name)        

        if run == 'training':
            self.io.print_info('Finished setting Training and Validation samples')
            self.n_collects = len(self.total_data)
            self.training_ids = list()
            self.val_ids = list()
            for i in range(self.n_collects):
                temp_n_sample = len(self.total_data[i])
                temp_ids_for_list = np.arange(temp_n_sample)
                np.random.shuffle(temp_ids_for_list)
                self.training_ids.append(temp_ids_for_list[:int(self.cfgs['train_split'] * temp_n_sample)])
                self.val_ids.append(temp_ids_for_list[int(self.cfgs['train_split'] * temp_n_sample):])
            self.tr_per_collect = self.cfgs[run]['trace_per_collect']
            self.rand_dev = self.cfgs[run]['rand_dev']
            self.rand_dev_ne = self.rand_dev * (-1)
            self.length_before = self.cfgs[run]['length_before']
        elif run == 'testing':
            pass
            
        else:
            self.io.print_info('ERROR: CONFUSED ABOUT RUNNING MODE')

    def get_training_batch(self):
        batch_ids = list()
        for i in range(len(self.training_ids)):
            batch_ids.append(np.random.choice(self.training_ids[i], self.cfgs['batch_size_train']))
        return self.get_batch(batch_ids, 'batch_size_train')

    def get_validation_batch(self):
        batch_ids = list()
        for i in range(len(self.val_ids)):
            batch_ids.append(np.random.choice(self.val_ids[i], self.cfgs['batch_size_val']))
        return self.get_batch(batch_ids, 'batch_size_val')
    
    def get_batch(self,batch_ids,data_type):                
        data_height = 1
        
        batch_size = len(batch_ids) * self.cfgs[data_type]
        
        if self.cfgs['training']['add_noise'] == 'Y':
            batch_size += 3
        
        waves = np.zeros([int(batch_size),
                           int(self.cfgs['training']['n_channels']),
                           int(data_height),
                           int(self.cfgs['training']['data_width'])])
        
        picks = np.zeros([int(batch_size),
                           1,
                           int(data_height),
                           int(self.cfgs['training']['data_width'])])
        idx = 0
        
        for mdx in range(len(batch_ids)):
            for ndx in range(self.cfgs[data_type]):
                temp_id = batch_ids[mdx][ndx]
                if self.cfgs['training']['n_channels'] > 1:
                    for cha_dx in range(self.cfgs['training']['n_channels']):
                        temp_wave = np.reshape(self.total_data[mdx][temp_id][0],[self.cfgs['training']['n_channels'],self.cfgs['training']['data_width']])
                        waves[idx,cha_dx,0,:] = temp_wave[cha_dx][:]
                else:                    
                    #print('{} of {}'.format(temp_id,mdx))
                    waves[idx,0,0,:] = self.total_data[mdx][temp_id][0][:]
                    
                temp_gt = np.zeros(self.cfgs['training']['data_width'])
                temp_gt[self.total_data[mdx][temp_id][1]] = 1.0
                temp_gt[self.total_data[mdx][temp_id][1] + 1] = 1.0
                temp_gt[self.total_data[mdx][temp_id][1] - 1] = 1.0
                picks[idx,0,0,:] = temp_gt[:]
                idx = idx + 1
        
        if self.cfgs['training']['add_noise'] == 'Y':
            waves[-1,0,0,:] = np.random.normal(0,0.5,self.cfgs['training']['data_width'])
            waves[-2,0,0,:],waves[-3,0,0,:] = self.get_random_triangle_noise()
            
        waves = np.reshape(waves,(int(batch_size),
                                   1,
                                   int(self.cfgs['training']['data_width']),
                                   int(self.cfgs['training']['n_channels']),
                                   ))
        
        picks = np.reshape(picks,(int(batch_size),
                                   1,
                                   int(self.cfgs['training']['data_width']),
                                   1,
                                   ))
        """
        np.save('wavescheck.npy',waves)
        np.save('pickscheck.npy',picks)
        """
        return waves, picks
    
    def get_random_triangle_noise(self):
        triangle_one_side = np.zeros(self.cfgs['training']['data_width'])
        triangle_two_side = np.zeros(self.cfgs['training']['data_width'])
        
        tri_num_one = np.random.randint(low = 10, high = 20)
        tri_num_two = np.random.randint(low = 10, high = 20)
        
        for i in range(tri_num_one):
            temp_gap = np.random.randint(30,50)
            temp_tri_width = np.random.randint(2,14)
            temp_gain = 0.5/temp_tri_width
            triangle_one_side[100+temp_gap*i+temp_tri_width] = 0.5
            for j in range(temp_tri_width):
                triangle_one_side[100+temp_gap*i+j] = temp_gain * j
                triangle_one_side[100+temp_gap*i+temp_tri_width*2 - j] = temp_gain * j
        flag = -1.0
        for i in range(tri_num_two):
            temp_gap = np.random.randint(30,50)
            temp_tri_width = np.random.randint(2,14)
            temp_gain = 0.7/temp_tri_width*flag
            triangle_one_side[100+temp_gap*i+temp_tri_width] = 0.7*flag
            flag = flag * (-1.0)
            for j in range(temp_tri_width):
                triangle_two_side[100+temp_gap*i+j] = temp_gain * j
                triangle_two_side[100+temp_gap*i+temp_tri_width*2 - j] = temp_gain * j
        
        return triangle_one_side,triangle_two_side