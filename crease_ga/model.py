import numpy as np
from os import path
import os
from crease_ga import utils
from crease_ga.adaptation_params import adaptation_params
import random    
import matplotlib
#matplotlib.use('Agg') ## uncomment this when running on cluster, comment out this line if on local
import matplotlib.pyplot as plt
import sys
from importlib import import_module
import time
from warnings import warn
from crease_ga.exceptions import CgaError

class Model:
    """
    The basic class that defines the model to be used to solve for a scattering
    profile.
    
    Attributes
    ----------
    pop_number: int. 
        Number of individuals within a generation.
    generations: int.
        Number of generations to run.
    nloci: int.
        Number of binary bits to represent each parameter in an individual.
        The decimal parameter value is converted to binary, with "all 0s"
        corresponding to the min value and "all 1s" corresponding to the
        max value. The larger the value, the finer the resultion for each parameter.
    adaptation_params: crease_ga.adaptation_params.adaptation_params
        Object of adaptation parameters used for this model.

    See also
    --------
    crease_ga.adaptaion_params.adaptation_params
    """

       
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
        '''
        Load a shape.

        Parameters
        ----------
        shape: str. name of the shape.
            Currently supported builtin shapes are "vesicle" and "micelle". Can
            also specify a shape developed in a crease_ga plugin.
        shape_params: list.
            Values of shape-specific descriptors. See the API of corresponding
            shape for details. If not specified, or an incorrect number of
            shape descriptor values are specified, the default values of the
            shape-specific descriptors will be loaded.
        minvalu,maxvalu: list.
            Values of the minimum and maximum boundaries of the
            parameters to be fit. If not specified, or an incorrect number of
            input parameter boundaries are specified, the default boundaries of
            the input parameters of the shape will be loaded.
        '''
        builtin_shapes=["vesicle","micelle","NP-solution","binary-NP-assembly",
                        "sffibril"]
        if shape in builtin_shapes:
            sg = import_module('crease_ga.shapes.'+shape+'.scatterer_generator')
            sg = sg.scatterer_generator
            print('imported builtin shape {}\n'.format(shape))
        else:
            from crease_ga.plugins import plugins
            if shape in plugins.keys():
                sg = plugins[shape].load()
                print('imported shape {} as a plugin'.format(shape))
            else:
                raise CgaError('Currently unsupported shape {}'.format(shape))
        
        #TODO: Complete the checker
        if shape_params is None:
            self.scatterer_generator = sg()
        elif minvalu == None or maxvalu == None:
            warn("Unspecified minimum and/or maximum parameter boundaries. Fall back to the default minimum "
                 "and maximum parameter boundaries of shape {}.\n".format(shape),stacklevel = 2)
            self.scatterer_generator = sg(shape_params)
            print("minimum parameter boundaries have been set to {},\n"
                  "maximum parameter boundaries have been set to {}.\n".format(
                   self.scatterer_generator.minvalu,
                   self.scatterer_generator.maxvalu))

        elif sg().numvars != len(minvalu) or sg().numvars != len(maxvalu):
               
            raise CgaError("Number of parameters in minvalu and/or maxvalu is not equal to number of parameters "
                 "required by shape {}.\n Shape {} requires {:d} parameters.\nminvalu has {:d} parameters.\n"
                 "maxvalu has {:d} parameters.".format(shape,shape,sg().numvars,
                                                     len(minvalu),len(maxvalu))) 
        else:
             self.scatterer_generator = sg(shape_params,minvalu,maxvalu)


        self.numvars = self.scatterer_generator.numvars   
        self.minvalu = self.scatterer_generator.minvalu
        self.maxvalu = self.scatterer_generator.maxvalu
            
            
            
    def load_iq(self,input_file_path,q_bounds=None):
        """
        Load an experimental I(q) profile [Iexp(q)] to the model, so that it can be
        solved later using "Model.solve".
        
        Parameters
        ----------
        input_file_path: str. Path to the input file. 
            The file should be organized in up to three column, with q-values in the first column, 
            corresponding I(q) values in the second, and optionally, corresponding error for the 
            I(q) in the third.
        q_bounds: [min,max].
            Define the minimum and maximum bound of the q region of interest. Any
            q-I(q) pairs outside of the defined bounds will be ignored during the
            fitting.
    
        See also
        --------
            crease_ga.Model.solve()
        """
        loadvals = np.genfromtxt(input_file_path)
        self.qrange_load = loadvals[:,0]
        IQin_load = loadvals[:,1]
        if len(loadvals.T)>2:
            IQerr_load = loadvals[:,2]
            IQerr_load = np.true_divide(IQerr_load,np.max(IQin_load))
        else:
            IQerr_load = None
        self.IQin_load=np.true_divide(IQin_load,np.max(IQin_load))
        #TODO: highQ and lowQ needs to be able to be dynamically set
        if q_bounds is None:
            self.qrange = self.qrange_load
            self.IQin = self.IQin_load
            self.IQerr = IQerr_load
        else:
            lowQ = q_bounds[0]
            highQ = q_bounds[1]
            self.IQin = self.IQin_load[ np.where(self.qrange_load>=lowQ)[0][0]:np.where(self.qrange_load<=highQ)[0][-1] +1]
            if IQerr_load is None:
                self.IQerr = None
            else:
                self.IQerr = IQerr_load[ np.where(self.qrange_load>=lowQ)[0][0]:np.where(self.qrange_load<=highQ)[0][-1] +1]
            self.qrange = self.qrange_load[ np.where(self.qrange_load>=lowQ)[0][0]:np.where(self.qrange_load<=highQ)[0][-1] +1]
            
        baseline = self.IQin[0]
        if self.IQerr is not None:
            self.IQerr = np.true_divide(self.IQerr,baseline)
        self.IQin = np.true_divide(self.IQin,baseline)

        
    def solve(self,name = 'ga_job',
              verbose = True,
              backend = 'debye',
              fitness_metric = 'log_sse',
              output_dir='./',
              n_cores=1,
              needs_postprocess = False):
        '''
        Fit the loaded target I(q) for a set of input parameters that maximize
        the fitness or minimize the error metric (fitness_metric).

        Parameters
        ----------
        name: str.
            Title of the current run. A folder of the name will be created
            under current working directory (output_dir), and all output files
            will be saved in that folder.
        verbose: bool. Default=True.
            If verbose is set to True, a figure will be produced at the end of
            each run, plotting the I(q) resulting from the best
            individual in the current generation and the target I(q).

            Useful for pedagogical purpose on jupyter notebook.
        fitness_metric: string. Default='log_sse'.
            The metric used to calculate fitness. Currently supported:
                "log_sse", sum of squared log10 difference at each q
                point.
        output_dir: string. Default="./" 
            Path to the working directory.
        '''
        ### checking if starting new run or restarting partial run
        address = output_dir+'/'+name+'/'
        if path.isfile(address+'current_gen.txt'):
            currentgen = int(np.genfromtxt(address+'current_gen.txt'))
            pop = np.genfromtxt(address+'current_pop.txt')
            temp = np.genfromtxt(address+'current_pm_pc.txt')
            pm = temp[0]
            pc = temp[1]
            self.adaptation_params.pc = pc
            self.adaptation_params.pm = pm
            # read in best iq for each generation
            bestIQ = np.genfromtxt(address+'best_iq.txt')
            # do not include q values in bestIQ array
            bestIQ = bestIQ[1:,:]
            print('Restarting from gen #{:d}'.format(currentgen+1))
        else:
            os.mkdir(address)
            currentgen = 0
            pop = utils.initial_pop(self.popnumber, self.nloci, self.numvars)
            # save best iq for each generation (plus q values)
            with open(address+'best_iq.txt','w') as f:
                np.savetxt(f,self.qrange,fmt="%-10f",newline='')
            bestIQ = []
            print('New run')
        
        colors = plt.cm.coolwarm(np.linspace(0,1,self.generations))
        for gen in range(currentgen, self.generations):    
            if backend == 'debye':
                pacc,gdm,elitei,IQid_str = self.fitness(pop,gen,output_dir+'/'+name+'/',fitness_metric,n_cores)
                IQid_str = np.array(IQid_str)
            pop = self.genetic_operations(pop,pacc,elitei)
            self.adaptation_params.update(gdm)
            if bestIQ == []:
                bestIQ = IQid_str[elitei]
            else:
                bestIQ = np.vstack((bestIQ,IQid_str[elitei]))
            with open(address+'best_iq.txt','a') as f:
                f.write('\n')
                np.savetxt(f,IQid_str[elitei],fmt="%-10f",newline='')

            ### save output from current generation in case want to restart run
            np.savetxt(address+'current_gen.txt',np.c_[gen])
            np.savetxt(address+'current_pop.txt',np.c_[pop])
            np.savetxt(address+'current_pm_pc.txt',np.c_[self.adaptation_params.pm,self.adaptation_params.pc])
            
            if needs_postprocess:
                self.postprocess()

            if verbose:
                figsize=(4,4)
                fig, ax = plt.subplots(figsize=(figsize))
                ax.plot(self.qrange_load,self.IQin,color='k',linestyle='-',ms=8,linewidth=1.3,marker='o')
                for i in range(gen+1):
                    ax.plot(self.qrange,bestIQ[i],color=colors[i],linestyle='-',ms=8,linewidth=2)
                plt.xlim(self.qrange[0],self.qrange[-1])
                plt.ylim(2*10**(-5),20)
                plt.xlabel(r'q, $\AA^{-1}$',fontsize=20)
                plt.ylabel(r'$I$(q)',fontsize=20)
                ax.set_xscale("log")
                ax.set_yscale("log")
                fig.savefig(output_dir+'plot'+str(gen)+'.png')
                plt.show()
                if gen == self.generations-1:
                    plt.savefig('iq_evolution.png',dpi=169,bbox_inches='tight')
    
    def postprocess(self):
        #import weakref
        self.scatterer_generator.postprocess(self)
      
    def fitness(self,pop,generation,output_dir,metric='log_sse',n_cores=1):
        tic = time.time()
        cs=10
        F1= open(output_dir+'results_'+str(generation)+'.txt','w')
        F1.write('#individual...all params...error\n')
        np.savetxt(output_dir+'population_'+str(generation)+'.txt',np.c_[pop])

        params=[]
        # calculate scattering for each individual
        for val in range(self.popnumber):
            param=utils.decode(pop, val, self.nloci, self.minvalu, self.maxvalu) # gets the current structure variables
            params.append(param)
        IQids = self.scatterer_generator.calculateScattering(self.qrange,params,output_dir,n_cores)
        
        fitn=np.zeros(self.popnumber)
        fitnfr=np.zeros(self.popnumber)
        fit=np.zeros(self.popnumber)
        qfin=self.qrange[-1]
        IQid_str=[]
        for val in range(self.popnumber):
            ### calculate computed Icomp(q) ###
            IQid=IQids[val]

            err=0
            for qi,qval in enumerate(self.qrange):
                if (IQid[qi]>0)&(self.IQin[qi]>0):
                    if (qi<qfin):
                        wil=np.log(np.true_divide(self.qrange[qi+1],self.qrange[qi]))  # weighting factor
                    else:
                        wil=np.log(np.true_divide(self.qrange[qi],self.qrange[qi-1]))  # weighting factor
                    if metric == 'log_sse':
                        err+=wil*(np.log(np.true_divide(self.IQin[qi],IQid[qi])))**2  # squared log error 
                elif metric == 'chi2':
                    if self.IQerr is None:
                        # chi^2 with weighting of IQin
                        err += np.true_divide(
                            np.square(self.IQin[qi]-IQid[qi]), np.square(self.IQin[qi]))
                    else:
                        # chi^2 with weighting of IQerr
                        err += np.true_divide(
                            np.square(self.IQin[qi]-IQid[qi]), np.square(self.IQerr[qi]))
            fit[val]=err
            IQid_str.append(IQid)  

            F1.write(str(val)+' ')
            for p in param:
                F1.write(str(p)+' ')
            F1.write(str(err)+'\n')
            F1.flush()
        self.fit = fit

        params = np.array(params)
        maxerr=np.max(fit)           #determines maximum SSerror for the population
        fitn=np.subtract(maxerr,fit) #determines error differences
        bestfit=np.max(fitn)
        sumup=np.sum(fitn)

        avgfit=np.true_divide(sumup,self.popnumber)
        dval=bestfit-avgfit
        ascale=np.true_divide(avgfit,dval)*(cs-1.0)     #linear scaling with cs as a scaleFactor
        bscale=avgfit*(1.0-ascale)

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
        print('Generation gdm: {:.3f}'.format(gdm))
        print('Generation best parameters '+str(params[elitei]))
        IQid_str = np.array(IQid_str)
        with open(output_dir+'IQid_best.txt','a') as f:
            f.write(np.array2string(IQid_str[elitei][0])+'\n')

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
        
        

        
