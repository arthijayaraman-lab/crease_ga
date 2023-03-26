import numpy as np
import numexpr as ne
from scipy import stats
#from crease_ga.exceptions import CgaError
from crease_ga import utils
from os import path
import multiprocessing as mp
from tensorflow import keras

## this is for slit height = 0 & slit width = value
def weightMatrix(qext,q,h):
    # using SasView as template 
    qedge = np.hstack([qext[0]-0.5*(qext[1]-qext[0]),0.5*(qext[1:]+qext[:-1]),qext[-1]+0.5*(qext[-1]-qext[-2]),])
    weight = np.zeros((len(q),len(qext)))
    for i, qi in enumerate(q):
        # assumes w = 0
        inx = 1. * ((qext >= qi-h) & (qext <= qi+h))
        absx = 1. * (qext < abs(qi-h)) if qi<h else 0.
        weight[i,:] = (inx + absx) * np.diff(qedge) / (2.*h)
    return weight.transpose()

## this is for slit height = value & slit width = 0
def weightMatrix2(qext,q,h):
    qedge = np.hstack([qext[0]-0.5*(qext[1]-qext[0]),0.5*(qext[1:]+qext[:-1]),qext[-1]+0.5*(qext[-1]-qext[-2]),])
    weight = np.zeros((len(q),len(qext)))
    for i, qi in enumerate(q):
        u_limit = np.sqrt(qi**2 + h**2)
        u_edges = qedge**2 - qi**2
        u_edges[qedge < abs(qi)] = 0.
        u_edges[qedge > u_limit] = u_limit**2 - qi**2
        weight[i,:] = np.diff(np.sqrt(u_edges))/h
    return weight.transpose()

def smearSlit(q,qext,iq,h=0.23969):
    #weight = weightMatrix(qext,q,h)
    weight = weightMatrix2(qext,q,h)
    iqsmear = np.dot(iq[None,:],weight).flatten()
    return iqsmear

class scatterer_generator:
    def __init__(self,
                 shape_params=[1, 1,20000,False,0,'','',''],
                 minvalu=(40,0.05,50,0.01,0.05,0.01,0,0,0,10),
                 maxvalu=(60,0.20,50,0.01,0.95,0.50,1,1,1,20)):
        
        # allow input of scattering differences between A and B type polymers
        sld1 = shape_params[0]
        sld2 = shape_params[1]
        self.sld1 = sld1
        self.sld2 = sld2
        # number of particles, sets box size
        N = shape_params[2]
        self.slit = shape_params[3]
        self.slit_h = shape_params[4]
        self.mlModelsq = shape_params[5]
        self.mlModelfq = shape_params[6]
        self.sqnorm = shape_params[7]
        self.smear = self.slit
        self.N = N  # Number of particles to use
        
        # genes are:
            # 1) average micelle diameter
            # 2) micelle lognormal dispersity
            # 3) solvent bead diameter -- for spacing micelles out
            # 4) solvent lognormal dispersity
            # 5) A:B core:shell split -- fraction of total micelle diameter
            # 6) micelle concentration (vol)
            # 7) gene related to aggregation of micelles
            # 8) gene ralted to sphericity of aggregate
            # 9) gene related to aggreagte spacing
            # 10) background scattering
        self._numvars = 10
        self.minvalu = minvalu
        self.maxvalu = maxvalu
            
        if path.isfile('current_sse.txt'):
            self.best_fit = np.genfromtxt('current_sse.txt')
        else:
            self.best_fit = np.inf

    @property
    def numvars(self):
        return self._numvars
    
    def calculateScattering(self,model,qrange,params,output_dir,n_cores=1):
        '''
        Determine each individual's computed scattering intensity profile.

        Parameters
        ----------
        qrange: numpy.array
            q values.
        params: numpy.array
            Decoded input parameters. See *Notes* section of the class
            documentation.
        n_cores: int
            number of CPU cores to use during multiprocessing

        Returns
        -------
        IQids: A numpy array holding each individual's I(q).
        '''
        self.qrange = qrange
        self.output_dir = output_dir
        self.params = params
        iqall = []
        sqall = []
        fqall = []
        pool = mp.Pool(n_cores)
        one = pool.starmap(self.doScattering, zip(
                range(len(params)), params))
        pool.close()
        pool.join()
        print('finished scattering')
        one = np.array(one,dtype=object)
        iqall = np.array(one[:, 0])
        sqall = one[:, 1]
        fqall = one[:, 2]
        self.fqall = fqall
        self.sqall = sqall

        return iqall
    
    def doScattering(self, individual, param):
        if self.smear == True:
            if self.slit == True:
                h = self.slit_h
                qmin, qmax = np.min(self.qrange-h), np.max(np.sqrt((self.qrange-h)**2))
                if qmin < 0:
                    qmin = self.qrange[0]*0.02
                log_delta_q = (len(self.qrange) - 1) / (np.log(self.qrange[-1]) - np.log(self.qrange[0]))
                nlow = int(np.ceil(log_delta_q * (np.log(self.qrange[0])-np.log(qmin))))
                qlow = np.logspace(np.log10(qmin), np.log10(self.qrange[0]), nlow+1)[:-1]
                nhigh = int(np.ceil(log_delta_q * (np.log(qmax)-np.log(self.qrange[-1]))))
                qhigh = np.logspace(np.log10(self.qrange[-1]), np.log10(qmax), nhigh+1)[1:]
                qext = np.concatenate([qlow,self.qrange,qhigh])

            # get form factor values
            model = keras.models.load_model(self.mlModelfq)
            xval = np.zeros(3)
            fq = np.zeros(len(qext))
            for ind, qval in enumerate(qext):
                xval[0] = np.log10(qval*param[0]/(2*np.pi))
                xval[1] = param[4]
                xval[2] = param[1]
                xinput = xval.reshape((1,len(xval)))
                # trained ANN to predict -log10(Icomp)
                output = model(xinput)
                temp = np.power(10,-output[:3])
                temp[0] *= self.sld1*self.sld1
                temp[1] *= self.sld2*self.sld2
                temp[2] *= self.sld1*self.sld2
                temp[-1] *= np.where(output[-1]<0,-1,1)
                fq[ind] = np.sum(temp)
            
            # get sq values
            model = keras.models.load_model(self.mlModelsq)
            xval = np.zeros(6)
            sq = np.zeros(len(qext))
            sqnorm = np.loadtxt(self.sqnorm)
            for ind, qval in enumerate(qext):
                xval[0] = np.log10(qval*param[0]/(2*np.pi))
                xval[1] = param[1]
                xval[2] = param[5]
                xval[3] = param[6]
                xval[4] = param[7]
                xval[5] = param[8]
                xinput = xval.reshape((1,len(xval)))
                # trained ANN to predict -log10(Icomp)
                sq[ind] = model(xinput)

            icomp = fq*fq*sq
            icompSmear = smearSlit(self.qrange, qext, icomp, h)
            #normalize Icomp to 1.0 at lowest q
            icompSmear = np.true_divide(icompSmear,icompSmear[0])
            # add background
            Background = 10**(-param[9])
            icompSmear += Background
            return icompSmear
        else:
            q = self.qrange
            # get form factor values
            model = keras.models.load_model(self.mlModelfq)
            xval = np.zeros(3)
            IQid = np.zeros((len(q),3))
            for ind, qval in enumerate(q):
                xval[0] = np.log10(qval*param[0]/(2*np.pi))
                xval[1] = param[4]
                xval[2] = param[1]
                xinput = xval.reshape((1,len(xval)))
                # trained ANN to predict -log10(Icomp)
                output = np.array(model(xinput))[0]
                IQid[ind,:] = np.power(10,-output)
          
            fqtot = IQid[:,0]
            fqcore = IQid[:,1]
            fqshell = IQid[:,2] 
            del IQid
    
            # get sq values
            model = keras.models.load_model(self.mlModelsq)
            xval = np.zeros(6)
            sq = np.zeros(len(q))
            # sqnorm contains the normalization values (average and std dev) for
            #   applying the sq ANN at the q of interest (instead of q trained on)
            sqnorm = np.loadtxt(self.sqnorm)
            sqave = np.interp(q*param[0]/(2*np.pi),sqnorm[:,0],sqnorm[:,1])
            sqstd = np.interp(q*param[0]/(2*np.pi),sqnorm[:,0],sqnorm[:,2])
            for ind, qval in enumerate(q):
                xval[0] = np.log10(qval*param[0]/(2*np.pi))
                xval[1] = param[1]
                xval[2] = param[5]
                xval[3] = param[6]
                xval[4] = param[7]
                xval[5] = param[8]
                xinput = xval.reshape((1,len(xval)))
                # trained ANN to predict -log10(Icomp)
                sq[ind] = model(xinput)*sqstd[ind] + sqave[ind]

            icomptot = np.abs(fqtot*fqtot*sq)
            icompcore = np.abs(fqcore*fqcore*sq)
            icompshell = np.abs(fqshell*fqshell*sq)
            icomptot = np.true_divide(icomptot,icomptot[0])
            icompcore = np.true_divide(icompcore,icompcore[0])
            icompshell = np.true_divide(icompshell,icompshell[0])
            # add background
            Background = 10**(-param[9])
            icomptot += Background
            icompcore += Background
            icompshell += Background
            icomp = np.append(icomptot,icompcore)
            icomp = np.append(icomp,icompshell)
            return icomp, fqtot, sq
        
    def postprocess(self,model):
        '''
        Save struture of best overall individual

        Parameters
        ----------
        model: object
            model.py self reference
        '''
        # save necessary information for structure
        if np.min(model.fit) < self.best_fit:
            self.best_fit = np.min(model.fit)
            val = np.where(model.fit == self.best_fit)[0][0]
            best_param = self.params[val].flatten()
            del self.params
            np.savetxt(self.output_dir+'current_fit.txt',np.c_[self.best_fit])
            np.savetxt(self.output_dir+'current_genes.txt',np.c_[best_param])
            # get diameter
            D1 = best_param[0]
            S1 = best_param[1]
            D2 = best_param[2]
            S2 = best_param[3]
            binmin1 = stats.lognorm.ppf(0.01, S1, scale=D1)
            binmax1 = stats.lognorm.ppf(0.99, S1, scale=D1)
            sizes1 = np.linspace(binmin1, binmax1, 11)
            binmin2 = stats.lognorm.ppf(0.01, S2, scale=D2)
            binmax2 = stats.lognorm.ppf(0.99, S2, scale=D2)
            sizes2 = np.linspace(binmin2, binmax2, 11)
            diameters = np.append(sizes1,sizes2)
            np.savetxt(self.output_dir+'current_diameters.txt',np.c_[diameters])
            np.savetxt(self.output_dir+'best_fq.txt',np.c_[self.fqall[val]])
            np.savetxt(self.output_dir+'best_sq.txt',np.c_[self.sqall[val]])
