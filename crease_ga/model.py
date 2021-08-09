import numpy as np
from os import path
from crease_ga import utils
from crease_ga.adaptation_params import adaptation_params
import random    
import matplotlib
#matplotlib.use('Agg') ## uncomment this when running on cluster, comment out this line if on local
import matplotlib.pyplot as plt
import sys
from importlib import import_module
import time
import sys
from warnings import warn
from crease_ga.exceptions import CgaError

class Model:
        
    def __init__(self,
                 pop_number = 5,
                 generations = 10,
                 nloci = 7,
                 yaml_file='x'):
        
        if path.isfile(yaml_file):
            pass
            #TODO: populate all input parameters with input from yaml files
        else:
            self.popnumber =  pop_number
            self.generations = generations
            self.nloci = nloci
            #TODO: check numvars is equal to length of minvalu and maxvalu
        self.adaptation_params = adaptation_params()  
    def load_shape(self,shape="vesicle", shape_params=None,minvalu=None,maxvalu=None):
        
        builtin_shapes=["vesicle"]
        if shape in builtin_shapes:
            sg = import_module('crease_ga.shapes.'+shape+'.scatterer_generator')
        else:
            raise CgaError('Currently unsupported shape {}'.format(shape))
        
        #TODO: Complete the checker
        if shape_params == None:
            self.scatterer_generator = sg.scatterer_generator()
        elif minvalu == None or maxvalu == None:
            warn("Unspecified minimum and/or maximum parameter boundaries. Fall back to the default minimum "
                 "and maximum parameter boundaries of shape {}.\n".format(shape),stacklevel = 2)
            self.scatterer_generator = sg.scatterer_generator(shape_params)
            print("minimum parameter boundaries have been set to {},\n"
                  "maximum parameter boundaries have been set to {}.\n".format(
                   self.scatterer_generator.minvalu,
                   self.scatterer_generator.maxvalu))

        elif sg.scatterer_generator().numvars != len(minvalu) or sg.scatterer_generator().numvars != len(maxvalu):
               
            raise CgaError("Number of parameters in minvalu and/or maxvalu is not equal to number of parameters "
                 "required by shape {}.\n Shape {} requires {:d} parameters.\nminvalu has {:d} parameters.\n"
                 "maxvalu has {:d} parameters.".format(shape,shape,sg.scatterer_generator().numvars,
                                                     len(minvalu),len(maxvalu))) 
        else:
             self.scatterer_generator = sg.scatterer_generator(shape_params,minvalu,maxvalu)


        self.numvars = self.scatterer_generator.numvars   
        self.minvalu = self.scatterer_generator.minvalu
        self.maxvalu = self.scatterer_generator.maxvalu
            
            
            
    def load_iq(self,input_file_path,q_bounds=None):
        loadvals = np.loadtxt(input_file_path)
        self.qrange_load = loadvals[:,0]
        IQin_load = loadvals[:,1]
        self.IQin_load=np.true_divide(IQin_load,np.max(IQin_load))
        #TODO: highQ and lowQ needs to be able to be dynamically set
        if q_bounds == None:
            self.qrange = self.qrange_load
            self.IQin = self.IQin_load
        else:
            lowQ = q_bounds[0]
            highQ = q_bounds[1]
            self.IQin = self.IQin_load[ np.where(self.qrange_load>=lowQ)[0][0]:np.where(self.qrange_load>=highQ)[0][0] ]
            self.qrange = self.qrange_load[ np.where(self.qrange_load>=lowQ)[0][0]:np.where(self.qrange_load>=highQ)[0][0] ]
            

        baseline = self.IQin[0]
        self.IQin = np.true_divide(self.IQin,baseline)
        self.IQin_load = np.true_divide(self.IQin_load,baseline)

        
    def solve(self,verbose = True,backend = 'debye',fitness_metric = 'log_sse',output_dir='./'):
        pop = utils.initial_pop(self.popnumber, self.nloci, self.numvars)
        for gen in range(self.generations):    
            if backend == 'debye':
                pacc,gdm,elitei,IQid_str = self.fitness(pop,gen,output_dir,metric='log_sse')
                IQid_str = np.array(IQid_str)
            pop = self.genetic_operations(pop,pacc,elitei)
            self.adaptation_params.update(gdm)
            
            if verbose:
                figsize=(4,4)
                fig, ax = plt.subplots(figsize=(figsize))
                ax.plot(self.qrange_load,self.IQin_load,color='k',linestyle='-',ms=8,linewidth=1.3,marker='o')
                ax.plot(self.qrange,IQid_str[elitei].transpose(),color='fuchsia',linestyle='-',ms=8,linewidth=2)#,marker='o')
                plt.xlim(0.001,0.1)
                plt.ylim(2*10**(-5),20)
                plt.xlabel(r'q, $\AA^{-1}$',fontsize=20)
                plt.ylabel(r'$I$(q)',fontsize=20)
                ax.set_xscale("log")
                ax.set_yscale("log")
                plt.show()

            
    def fitness(self,pop,generation,output_dir,metric='log_sse'):
        tic = time.time()
        cs=10
        F1= open(output_dir+'results_'+str(generation)+'.txt','w')
        F1.write('#individual...all params...error\n')
        np.savetxt(output_dir+'population_'+str(generation)+'.txt',np.c_[pop])
        
        fitn=np.zeros(self.popnumber)
        fitnfr=np.zeros(self.popnumber)
        fit=np.zeros(self.popnumber)
        qfin=self.qrange[-1]
        IQid_str=[]
        params=[]
        for val in range(self.popnumber):
            sys.stdout.write("\rGen {:d}/{:d}, individual {:d}/{:d}\n".format(generation+1,self.generations,val+1,self.popnumber))
            sys.stdout.flush()
            param=utils.decode(pop, val, self.nloci, self.minvalu, self.maxvalu) # gets the current structure variables
            params.append(param)

            ### calculate computed Icomp(q) ###
            IQid=self.scatterer_generator.converttoIQ(self.qrange, param)

            err=0
            for qi,qval in enumerate(self.qrange):
                if (IQid[qi]>0)&(self.IQin[qi]>0):
                    if (qi<qfin):
                        wil=np.log(np.true_divide(self.qrange[qi+1],self.qrange[qi]))  # weighting factor
                    else:
                        wil=np.log(np.true_divide(self.qrange[qi],self.qrange[qi-1]))  # weighting factor
                    if metric == 'log_sse':
                        err+=wil*(np.log(np.true_divide(self.IQin[qi],IQid[qi])))**2  # squared log error 
            fit[val]=err
            IQid_str.append(IQid)

            

            F1.write(str(val)+' ')
            for p in param:
                F1.write(str(p)+' ')
            F1.write(str(err)+'\n')
            F1.flush()

        params = np.array(params)
        maxerr=np.max(fit)           #determines maximum SSerror for the population
        fitn=fit-maxerr #determines error differences
        #maxfit=np.max(fitn)
        sumup=np.sum(fitn)

        avgfit=np.true_divide(sumup,self.popnumber)
        dval=-avgfit
        ascale=np.true_divide(avgfit,dval)*(cs-1.0)     #linear scaling with cs as a scaleFactor
        bscale=avgfit*(1.0-ascale)

        sumup=0

        # get scaled fitness to enable selection of bad candidates
        for val in range(self.popnumber):
            if (fitn[val]>avgfit):
                fitnfr[val]=ascale*fitn[val]+bscale
            else:
                fitnfr[val]=fitn[val]

        sumup=np.sum(fitnfr)

        pacc=np.zeros(self.popnumber)
        prob=np.true_divide(fitnfr,sumup)
        pacc=np.cumsum(prob)

        ### returns cummulative relative error from which individuals can be selected ###
        maxfit=np.min(fit)
        elitei=np.where(fit==maxfit)[0]                  # Best candidate 
        secondfit=sorted(fit)[1]
        secondi = np.where(fit==secondfit)[0]            # Second best candidate
        avgfit=np.average(fit)
        avgi=np.array([(np.abs(fit-avgfit)).argmin()])   # Average candidate
        minfit=np.max(fit)
        mini=np.where(fit==minfit)[0]                    # Worst candidate
        if avgfit==0:
            avgfit=1
        gdm=np.true_divide(maxfit,avgfit)
        if len(elitei)>1:
            elitei=elitei[0]
        if len(secondi)>1:
            secondi=secondi[0]
        if len(avgi)>1:
            avgi=avgi[0]
        if len(mini)>1:
            mini=mini[0]
        
        f = open(output_dir+'fitness_vs_gen.txt', 'a' )
        if generation == 0:
            f.write( 'gen mini min avgi avg secondi second besti best\n' )
        f.write( '%d ' %(generation) )
        f.write( '%d %.8lf ' %(mini,minfit) )
        f.write( '%d %.8lf ' %(avgi,avgfit) )
        f.write( '%d %.8lf ' %(secondi,secondfit) )
        f.write( '%d %.8lf ' %(elitei,maxfit) )
        f.write( '\n' )
        f.close()
        print('Generation time: {:.3f}s'.format(time.time()-tic))
        print('Generation best fitness: {:.4f}'.format(maxfit))
        print('Generation best fitness: {:.3f}'.format(gdm))
        print('Generation best parameters '+str(params[elitei]))
        
        return pacc, gdm, elitei, IQid_str
        
    
    def genetic_operations(self,pop,pacc,elitei):
        popn = np.zeros(np.shape(pop))
        cross = 0
        mute = 0
        pc = self.adaptation_params.pc
        pm = self.adaptation_params.pm
        
        for i in range(self.popnumber-1):

            #####################    Crossover    ####################

            #Selection based on fitness
            testoff=random.random()
            isit=0
            npart1=1
            for j in range(1,self.popnumber):
                if (testoff>pacc[j-1])&(testoff<pacc[j]):
                    npart1=j

            testoff=random.random()
            isit=0
            npart2=1
            for j in range(self.popnumber):
                if (testoff>=pacc[j-1])&(testoff!=pacc[j]):
                    npart2=j

            #Fit parents put in array popn

            popn[i,:]=pop[npart1,:]


            testoff=random.random()
            loc=int((testoff*(self.numvars-1))*self.nloci)
            if loc==0:
                loc=self.nloci
            testoff=random.random()

            #crossover
            if (testoff<=pc):
                cross+=1
                popn[i,loc:]=pop[npart2,loc:]


        #####################    Mutation    ####################


            for j in range(self.nloci*self.numvars):
                testoff=random.random()
                if (testoff<=pm):
                    popn[i,j]=random.randint(0,1)
                    mute+=1


        #####################    Elitism    ####################

        popn[-1,:]=pop[elitei,:]

        
        print('pc',pc)
        print('#crossovers',cross)
        print('pm',pm)
        print('#mutations',mute)
        print('\n')
        
        return popn
        
        

        
