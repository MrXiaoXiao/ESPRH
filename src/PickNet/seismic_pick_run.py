# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 08:58:55 2018
"""
import os
import argparse
import tensorflow as tf

from fcn.trainer import FCNTrainer
from fcn.tester import FCNTester

def get_session(gpu_fraction = None, gpuid = 0):
    num_threads = os.environ.get('OMP_NUM_THREADS')
    if num_threads is None:
        return tf.Session(config=tf.ConfigProto())
    else:
        num_threads = int(num_threads)
        
    gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=gpu_fraction)

    if num_threads:
        return tf.Session(config=tf.ConfigProto(gpu_options=gpu_options, intra_op_parallelism_threads=num_threads))
    else:
        return tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))

def main(args):
    os.environ['CUDA_VISIBLE_DEVICES'] = '{}'.format(args.gpuid)
    
    if not (args.run_train or args.run_test):
        print ('Please set one of the options --train | --test')
        parser.print_help()
        return
    
    if args.run_test or args.run_train:
        print('Setting GPU Limit')
        session = get_session(args.gpu_limit,args.gpuid)

    if args.run_train:
        print('Setting Trainer')
        trainer = FCNTrainer(args.config_file)
        trainer.setup()
        trainer.run(session)
        print('Done Training')

    if args.run_test:
        print('Setting Tester')
        tester = FCNTester(args.config_file)
        tester.setup(session)
        tester.run(session)
        print('Done Testing')
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Utility for Training/Testing PickNet models with Seismic Data')
    parser.add_argument('--config-file', dest='config_file', type=str, help='Experiment configuration file')
    parser.add_argument('--train', dest='run_train', action='store_true', default=False, help='Launch training')
    parser.add_argument('--test', dest='run_test', action='store_true', default=False, help='Launch testing on a list of images')
    parser.add_argument('--gpuid', dest='gpuid', type=int, default=0, help='Run On a certain GPU')
    parser.add_argument('--gpu-limit', dest='gpu_limit', type=float, default=1.0, help='Use fraction of GPU memory (Useful with TensorFlow backend)')
    args = parser.parse_args()
    main(args)        
