# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:11:05 2019

@author: michi
"""
'''
notes
popnumber = number of individuals
pop stores pop and genes
nloci = 7* #variables for the binary representation
'''
import numpy as np
import time
import multiprocessing as mp
from functools import partial
import sys
from crease_ga.shapes.ml_sffibril.nn_util import *
from tensorflow import keras
import tensorflow
from crease_ga.shapes.ml_sffibril.schulz_rv import *
from os import path
class scatterer_generator:
    '''
    The wrapper class for ML-enchanced semiflexible fibril shape. Default length unit: Angstrom.
    

    Notes
    -----
    **The following shape-specific descriptors are to be specified by user (see
    *Attributes*) as 
    a list, in the precise order as listed, while calling `Model.load_shape`
    to load this shape:**


    dist: 
        takes "constant" or "proportional".
        "constant": all monodisperse fibrils contain same number of
        "scatterers" (equivalent to same number of chains or monomers)
        regardless of their diameters. Thicker fibrils will have lower density.
        
        "proportional": all monodisperse fibrils have same packing density of
        "scatterers". Thicker chains will have more scatterers (more
        monomers or chains)
    sample_size:
        how many monodisperse fibril sampels are used to represent the
        polydispersity in diameter.
    model_path:
        path to the tensorflow model.
    db_ranges:
        min values,max values and qranges for normalization of inputs for the tensorflow model
        



    **The following 5 parameters are to be predicted, in the precise order
    as listed, by GA:**
    l_mu: mean fibril length. default [min,max]:[1000,4000]
    d_mu: mean fibril diameter. default [min,max]:[100,250]
    sigma_mu: mean fibril flexibility factor. default [min,max]: [0.05,0.8]
    d_s: standard deviation of fibril diameter. default [min,max]: [0.02,0.5]
    [min,max]:[0,1e-5]
    bg_log10: log10 of background scattering. default [min,max]:[-6,-2]
    '''
    def __init__(self,
        shape_params = ['constant',100,path.dirname(__file__)+'/sffibril_model.h5',
            [[200,50,0.005,-6],[5000,300,0.8,-2],[0.00324004,0.1346058]]],
        minvalu = [1000,100,0.05,0.1,-6],
        maxvalu = [4000,250,0.80,0.5,-2]):
    
        self._numvars = 5
    
        
        self.minvalu = minvalu
        self.maxvalu = maxvalu
        self.dist = shape_params[0]
        self.sample_size = shape_params[1]
        self.db = nn_db(shape_params[3][2],shape_params[3][0],shape_params[3][1])
        self.db.model = keras.models.load_model(shape_params[2])
    @property
    def numvars(self):
        return self._numvars
    
    def gen_eval_iq(self,val,params,qrange):
        param = params[val]
        l_mu = param[0]
        d_mu = param[1]
        sigma_mu = param[2]
        d_s = param[3]
        bg_log10 = param[4]
        ind_params = []
        weights = []
        IQids = []
        delta_d = (400-50)/self.sample_size
        for i in range(self.sample_size):
            
            di = delta_d*i+50
            l = l_mu
            s = sigma_mu
            ind_params.append(np.array([l,di,s,bg_log10]))
            weights.append(schulz_indpoint(di,d_mu,d_s,eval_method='log'))
        for ind_param in ind_params:

            xs = np.array([np.append(ind_param,q) for q in qrange])
            xs = self.db.preprocess(xs)
            m = self.db.model
            if self.dist == 'constant':
                IQids.append(np.array([10**i for i in m(xs).numpy()]))
            if self.dist == 'proportional':
                IQids.append(np.array([10**i for i in m(xs).numpy()])*ind_param[1]**4)
        weights = np.array(weights)
        ret = np.average(IQids,axis=0,weights=weights/np.sum(weights))
        return ret.flatten()


    def calculateScattering(self,qrange,params,output_dir,n_cores=1):
        '''
        Returns
        -------
        IQids: A numpy array holding each individual's I(q).
        '''
        if path.isfile(output_dir+'current_sse.txt'):
            self.best_fit = np.genfromtxt(output_dir+'current_sse.txt')
        else:
            self.best_fit = 1e10
        ret = [self.gen_eval_iq(val,params,qrange) for val in range(len(params))]
#        pool = mp.Pool(n_cores)
#        partial_work = partial(self.gen_eval_iq,
#                               params = params,
#                               qrange = qrange,
#                              )
#        IQids = pool.map(partial_work,[val for val in range(len(params))])
#        pool.close()
#        pool.join()
        return ret

    def postprocess(self,model):
        pass
