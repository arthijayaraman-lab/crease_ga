import numpy as np
import random
import numexpr as ne
from crease_ga.exceptions import CgaError
import sys
from crease_ga import utils
from subprocess import check_output, run
from itertools import repeat
from os import path
import multiprocessing as mp
from crease_ga.shapes.sffibril.sf_util import *
from functools import partial


class scatterer_generator:
    '''
    The wrapper class for semiflexible fibril shape. Default length unit: Angstrom.
    

    Notes
    -----
    **The following shape-specific descriptors are to be specified by user (see
    *Attributes*) as 
    a list, in the precise order as listed, while calling `Model.load_shape`
    to load this shape:**


    scatterer_density: 
        Density of scatterers used to represent the fibril, in A^-3. Default:
        0.5
    L:
        Default length of fibril. Default: 2000 A



    **The following 4 parameters are to be predicted, in the precise order
    as listed, by GA:**
    
    D_mean:
        mean fibril diameter. Default [min,max]: [60 A, 300 A]
    PD_D:
        Standard_Deviation_in_D/D_mean. Default [min,max]: [0.02, 0.5]
    sigma_theta:
        parameter describing chain stiffness. See paper for detail.
        Default [min,max]: [0.005, 0.6]
    -log10(bg):
        Negative log10 of background intensity. 
        E.g. an background intensity of 0.001 leads to this value being 3.
        Default [min,max]:  [2,6]
        '''
    def __init__(self,
                 shape_params=[2000, 0.0005],
                 minvalu=(60,0.02,0.005,2),
                 maxvalu=(300,0.5,0.6,6), 
                 custom_form=None):
        self._numvars = 4
        self.L = shape_params[0]
        self.rho_scat = shape_params[1]
        self.minvalu = minvalu
        self.maxvalu = maxvalu

    @property
    def numvars(self):
        return self._numvars
    
    def gen_eval_iq(self,val,params,qrange,sigmabead):
        #print("val="+str(val))
        param = params[val]
        dm = param[0]
        ds = param[1]*dm
        sample_size = 10
        IQid = np.zeros(len(qrange))
        d_array = np.linspace(dm-ds*1.5,dm+ds*1.5,sample_size)
        for d in d_array:
            scat_coords = gen_scat_coords_flexcyl(int(self.rho_scat*self.L*3.1416*(d/2)**2),2000,d,param[2])
            IQid_curr = sq(scat_coords,qrange)*np.array([spherical_pq(sigmabead,q) for q in qrange])
            normal_factor = 1/(ds*np.sqrt(2*np.pi))*np.exp(-(d-dm)**2/(2*ds**2))
            IQid += IQid_curr*normal_factor/sample_size # plus background

        return IQid/IQid[0]+10**(-param[3])    
    
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
        pool = mp.Pool(n_cores)
        partial_work = partial(self.gen_eval_iq,
                               params = params,
                               qrange = qrange,
                               sigmabead = 0.1)
        IQids = pool.map(partial_work,[val for val in range(len(params))])
        pool.close()
        pool.join()
        return IQids

    def postprocess(self,model):
        pass
