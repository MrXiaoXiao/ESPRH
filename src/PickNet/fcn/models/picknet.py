# Adapted from : VGG 16 model : https://github.com/machrisaa/tensorflow-vgg
import time
import os
import inspect

import numpy as np
from termcolor import colored
import tensorflow as tf

from fcn.losses import sigmoid_cross_entropy_balanced
from fcn.utils.io import IO


class PickNet():

    def __init__(self, cfgs, run='training'):

        self.cfgs = cfgs
        self.io = IO()

        self.inputs = tf.placeholder(tf.float32, [None, None, self.cfgs[run]['data_width'], self.cfgs[run]['n_channels']])
        self.groundtruths = tf.placeholder(tf.float32, [None , None, self.cfgs[run]['data_width'], 1])

        self.define_model()

    def define_model(self):

        """
        PickNet model based on RSRN
        Add branch layers (with deconv) after each CONV block
        """

        start_time = time.time()
        
         # block 1
        self.conv1_1 = self.conv_layer_vgg(self.inputs, name="conv1_1", 
                                           kh=self.cfgs['b1_convh'], kw=self.cfgs['b1_convw'], 
                                           n_out=64, dh=1, dw=1)
        self.conv1_2 = self.conv_layer_vgg(self.conv1_1,  name="conv1_2",
                                           kh=self.cfgs['b1_convh'], kw=self.cfgs['b1_convw'],  
                                           n_out=64, dh=1, dw=1)
        
        self.pool1 = self.max_pool(self.conv1_2,   name="pool1",   kh=1, kw=2, dh=1, dw=2)
        
        # block 2
        self.conv2_1 = self.conv_layer_vgg(self.pool1,    name="conv2_1", 
                                           kh=self.cfgs['b2_convh'], kw=self.cfgs['b2_convw'], 
                                           n_out=128, dh=1, dw=1)
        self.conv2_2 = self.conv_layer_vgg(self.conv2_1,  name="conv2_2",
                                           kh=self.cfgs['b2_convh'], kw=self.cfgs['b2_convw'],  
                                           n_out=128, dh=1, dw=1)
        
        self.pool2 = self.max_pool(self.conv2_2,   name="pool2",   kh=1, kw=2, dh=1, dw=2)
        
        # # block 3
        self.conv3_1 = self.conv_layer_vgg(self.pool2,    name="conv3_1", 
                                           kh=self.cfgs['b3_convh'], kw=self.cfgs['b3_convw'], 
                                           n_out=256, dh=1, dw=1)
        self.conv3_2 = self.conv_layer_vgg(self.conv3_1,  name="conv3_2", 
                                           kh=self.cfgs['b3_convh'], kw=self.cfgs['b3_convw'], 
                                           n_out=256, dh=1, dw=1)
        self.conv3_3 = self.conv_layer_vgg(self.conv3_2,  name="conv3_3", 
                                           kh=self.cfgs['b3_convh'], kw=self.cfgs['b3_convw'], 
                                           n_out=256, dh=1, dw=1)    
       
        self.pool3 = self.max_pool(self.conv3_3,   name="pool3",   kh=1, kw=2, dh=1, dw=2)
        
        # block 4
        self.conv4_1 = self.conv_layer_vgg(self.pool3,    name="conv4_1", 
                                           kh=self.cfgs['b4_convh'], kw=self.cfgs['b4_convw'], 
                                           n_out=512, dh=1, dw=1)
        self.conv4_2 = self.conv_layer_vgg(self.conv4_1,  name="conv4_2", 
                                           kh=self.cfgs['b4_convh'], kw=self.cfgs['b4_convw'], 
                                           n_out=512, dh=1, dw=1)
        self.conv4_3 = self.conv_layer_vgg(self.conv4_2,  name="conv4_3", 
                                           kh=self.cfgs['b4_convh'], kw=self.cfgs['b4_convw'],  
                                           n_out=512, dh=1, dw=1)
        
        self.pool4 = self.max_pool(self.conv4_3,   name="pool4",   kh=1, kw=2, dh=1, dw=2)
        
        # block 5
        self.conv5_1 = self.conv_layer_vgg(self.pool4,    name="conv5_1", 
                                           kh=self.cfgs['b5_convh'], kw=self.cfgs['b5_convw'],  
                                           n_out=512, dh=1, dw=1)
        self.conv5_2 = self.conv_layer_vgg(self.conv5_1,  name="conv5_2", 
                                           kh=self.cfgs['b5_convh'], kw=self.cfgs['b5_convw'],
                                           n_out=512, dh=1, dw=1)
        self.conv5_3 = self.conv_layer_vgg(self.conv5_2,  name="conv5_3", 
                                           kh=self.cfgs['b5_convh'], kw=self.cfgs['b5_convw'], 
                                           n_out=512, dh=1, dw=1)
        
        self.side_1_1 = self.side_layer(self.conv1_1,'side_1_1')
        self.side_1_2 = self.side_layer(self.conv1_2,'side_1_2')
        
        self.side_2_1 = self.side_layer(self.conv2_1,'side_2_1')
        self.side_2_2 = self.side_layer(self.conv2_2,'side_2_2')
        
        self.side_3_1 = self.side_layer(self.conv3_1,'side_3_1')
        self.side_3_2 = self.side_layer(self.conv3_2,'side_3_2')
        self.side_3_3 = self.side_layer(self.conv3_3,'side_3_3')

        self.side_4_1 = self.side_layer(self.conv4_1,'side_4_1')
        self.side_4_2 = self.side_layer(self.conv4_2,'side_4_2')
        self.side_4_3 = self.side_layer(self.conv4_3,'side_4_3')        
        
        self.side_5_1 = self.side_layer(self.conv5_1,'side_5_1')
        self.side_5_2 = self.side_layer(self.conv5_2,'side_5_2')
        self.side_5_3 = self.side_layer(self.conv5_3,'side_5_3')  
        
        #RSRN in a deep to shallow fashion
        #block 5
        self.dsnout_5_3 = self.deconv_layer(x = self.side_5_3,upscale=16,
                                            name='{}_dsnout_{}'.format('side_5_3', 16),
                                            w_init=tf.truncated_normal_initializer(stddev=0.1))
        
        self.residual_5_2, self.dsnout_5_2 = self.RU_layer(side_inputs = self.side_5_2, 
                                                           res_inputs = self.side_5_3, 
                                                           name = 'RU_5_2', 
                                                           upscale_out = 16, 
                                                           upscale_res = 1)
        
        self.residual_5_1, self.dsnout_5_1 = self.RU_layer(self.side_5_1, self.residual_5_2,'RU_5_1', 16, 1)
        
        #block 4
        self.residual_4_3, self.dsnout_4_3 = self.RU_layer(self.side_4_3, self.residual_5_1,'RU_4_3', 8, 2)
        self.residual_4_2, self.dsnout_4_2 = self.RU_layer(self.side_4_2, self.residual_4_3,'RU_4_2', 8, 1)
        self.residual_4_1, self.dsnout_4_1 = self.RU_layer(self.side_4_1, self.residual_4_2,'RU_4_1', 8, 1)
        
        #block 3
        self.residual_3_3, self.dsnout_3_3 = self.RU_layer(self.side_3_3, self.residual_4_1,'RU_3_3', 4, 2)
        self.residual_3_2, self.dsnout_3_2 = self.RU_layer(self.side_3_2, self.residual_3_3,'RU_3_2', 4, 1)
        self.residual_3_1, self.dsnout_3_1 = self.RU_layer(self.side_3_1, self.residual_3_2,'RU_3_1', 4, 1)
        
        #block 2
        self.residual_2_2, self.dsnout_2_2 = self.RU_layer(self.side_2_2, self.residual_3_1,'RU_2_2', 2, 2)
        self.residual_2_1, self.dsnout_2_1 = self.RU_layer(self.side_2_1, self.residual_2_2,'RU_2_1', 2, 1)
        
        #block 1
        self.residual_1_2, self.dsnout_1_2 = self.RU_layer(self.side_1_2, self.residual_2_1,'RU_1_2', 1, 2)
        self.residual_1_1, self.dsnout_1_1 = self.RU_layer(self.side_1_1, self.residual_1_2,'RU_2_1', 1, 1)
                   
        self.side_outputs = [self.dsnout_1_1, self.dsnout_1_2,
                             self.dsnout_2_1, self.dsnout_2_2,
                             self.dsnout_3_1, self.dsnout_3_2, self.dsnout_3_3,
                             self.dsnout_4_1, self.dsnout_4_2, self.dsnout_4_3,
                             self.dsnout_5_1,self.dsnout_5_2,self.dsnout_5_3]
        """
        self.side_outputs = [self.dsnout_1_1, self.dsnout_1_2]
        """ 
        
        w_shape = [1, 1, len(self.side_outputs), 1]
        
        self.fuse = self.conv_layer(tf.concat(self.side_outputs, axis=3),
                                    w_shape, name='fuse_1', use_bias=False,
                                    w_init=tf.constant_initializer(0.2))    
        
        
        self.outputs = self.side_outputs + [self.fuse]
        
        self.io.print_info("Build model finished: {:.4f}s".format(time.time() - start_time))
            
    
    def max_pool(self, input_tenosr, name,kh,kw,dh,dw):
        return tf.nn.max_pool(input_tenosr,
                              ksize=[1, kh, kw, 1],
                              strides=[1, dh, dw, 1], 
                              padding='SAME',name=name)
    
    def conv_layer_vgg(self,input_tensor, name,kh,kw,n_out,dh,dw):
        n_in = input_tensor.get_shape().as_list()[3]
        with tf.variable_scope(name):
            filt = tf.get_variable("weight",[kh,kw,n_in,n_out],
                                   initializer=tf.truncated_normal_initializer(stddev=0.1))
            conv = tf.nn.conv2d(input_tensor, filt, [1, dh, dw, 1], padding='SAME')
            conv_biases = tf.get_variable("bias",[n_out],initializer=tf.constant_initializer(0.0))
            bias = tf.nn.bias_add(conv, conv_biases)
            relu = tf.nn.relu(bias)
            return relu
    

    def conv_layer(self, x, W_shape, b_shape=None, name=None,
                   padding='SAME', use_bias=False, w_init=None, b_init=None):

        W = self.weight_variable(W_shape, w_init)
        tf.summary.histogram('weights_{}'.format(name), W)

        if use_bias:
            b = self.bias_variable([b_shape], b_init)
            tf.summary.histogram('biases_{}'.format(name), b)

        conv = tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding=padding)

        return conv + b if use_bias else conv

    def deconv_layer(self, x, upscale, name, padding='SAME', w_init=None):

        x_shape = tf.shape(x)
        in_shape = x.shape.as_list()

        w_shape = [1, upscale * 2, in_shape[-1], 1]
        strides = [1, 1, upscale, 1]

        W = self.weight_variable(w_shape, w_init)
        tf.summary.histogram('weights_{}'.format(name), W)

        out_shape = tf.stack([x_shape[0], x_shape[1], x_shape[2], w_shape[2]]) * tf.constant(strides, tf.int32)
        deconv = tf.nn.conv2d_transpose(x, W, out_shape, strides=strides, padding=padding)

        return deconv
    
    def side_layer(self, inputs, name):
        with tf.variable_scope(name):
            in_shape = inputs.shape.as_list()
            w_shape = [1, 1, in_shape[-1], 1]
            
        classifier = self.conv_layer(inputs, w_shape, b_shape=1,
                             w_init=tf.constant_initializer(),
                             b_init=tf.constant_initializer(),
                             name=name + '_side')
        return classifier
    
    def RU_layer(self, side_inputs, res_inputs, name, upscale_out, upscale_res):
        res_inputs_deconv = self.deconv_layer(x = res_inputs,
                                              upscale=upscale_res,
                                              name='{}_RUdeconv_{}'.format(name, upscale_res),
                                              w_init=tf.truncated_normal_initializer(stddev=0.1))
        
        
        classifier = tf.concat(values=[side_inputs,res_inputs_deconv], axis=3, name = name + '_concat')
        
        
        with tf.variable_scope(name):
            #in_shape = classifier.shape.as_list()
            w_shape = [1, 3, 2, 1]
        
        res_outputs = self.conv_layer(x = classifier, W_shape = w_shape, b_shape=1,
                                      w_init=tf.constant_initializer(),
                                      b_init=tf.constant_initializer(),
                                      name=name + '_RU')
        side_outputs = self.deconv_layer(x = res_outputs,
                                         upscale=upscale_out,
                                         name='{}_dsnout_{}'.format(name, upscale_out),
                                         w_init=tf.truncated_normal_initializer(stddev=0.1))
        
        return res_outputs, side_outputs
    """
    def side_layer(self, inputs, name, upscale, slice_num):
        with tf.variable_scope(name):
            in_shape = inputs.shape.as_list()
            w_shape = [1, 1, in_shape[-1], slice_num]
            
        classifier = self.conv_layer(inputs, w_shape, b_shape=1,
                                     w_init=tf.constant_initializer(),
                                     b_init=tf.constant_initializer(),
                                     name=name + '_reduction')
        
        classifier = self.deconv_layer(classifier,
                                       upscale=upscale,
                                       name='{}_deconv_{}'.format(name, upscale),
                                       w_init=tf.truncated_normal_initializer(stddev=0.1))
    
        return classifier
    """
    def get_conv_filter(self, name):
        return tf.constant(self.data_dict[name][0], name="filter")

    def get_bias(self, name):
        return tf.constant(self.data_dict[name][1], name="biases")

    def weight_variable(self, shape, initial):

        init = initial(shape)
        return tf.Variable(init)

    def bias_variable(self, shape, initial):

        init = initial(shape)
        return tf.Variable(init)

    def setup_testing(self, session):
        self.predictions = []

        for idx, b in enumerate(self.outputs):
            output = tf.nn.sigmoid(b, name='output_{}'.format(idx))
            #output = tf.nn.relu(b, name='output_{}'.format(idx))
            self.predictions.append(output)

    def setup_training(self, session):

        """
            Apply sigmoid non-linearity to side layer ouputs + fuse layer outputs
            Compute total loss := side_layer_loss + fuse_layer_loss
            Compute predicted edge maps from fuse layer as pseudo performance metric to track
        """

        self.predictions = []
        self.loss = 0

        self.io.print_warning('Deep supervision application set to {}'.format(self.cfgs['deep_supervision']))

        for idx, b in enumerate(self.side_outputs):
            output = tf.nn.sigmoid(b, name='output_{}'.format(idx))
            cost = sigmoid_cross_entropy_balanced(b, self.groundtruths, name='cross_entropy{}'.format(idx))

            self.predictions.append(output)
            if self.cfgs['deep_supervision']:
                self.loss += (self.cfgs['loss_weights'] * cost)

        fuse_output = tf.nn.sigmoid(self.fuse, name='fuse')
        fuse_cost = sigmoid_cross_entropy_balanced(self.fuse, self.groundtruths, name='cross_entropy_fuse')

        self.predictions.append(fuse_output)
        self.loss += (self.cfgs['loss_weights'] * fuse_cost)

        pred = tf.cast(tf.greater(fuse_output, 0.5), tf.int32, name='predictions')
        error = tf.cast(tf.not_equal(pred, tf.cast(self.groundtruths, tf.int32)), tf.float32)
        self.error = tf.reduce_mean(error, name='pixel_error')

        tf.summary.scalar('loss', self.loss)
        tf.summary.scalar('error', self.error)

        self.merged_summary = tf.summary.merge_all()

        self.train_writer = tf.summary.FileWriter(self.cfgs['save_dir'] + '/train', session.graph)
        self.val_writer = tf.summary.FileWriter(self.cfgs['save_dir'] + '/val')