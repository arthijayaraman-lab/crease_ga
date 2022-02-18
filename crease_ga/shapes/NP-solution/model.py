import numpy as np
from os import path
from crease_ga_diameter import utils
from crease_ga_diameter.adaptation_params import adaptation_params
import random
#import matplotlib
# matplotlib.use('Agg') ## uncomment this when running on cluster, comment out this line if on local
import matplotlib.pyplot as plt
import sys
from importlib import import_module
import time
from warnings import warn
from crease_ga_diameter.exceptions import CgaError
import multiprocessing as mp
from subprocess import run, check_output
from itertools import repeat

mp.process.current_process()._config['tempdir'] = '$(TMPDIR}'

class Model:

    def __init__(self,
                 pop_number=5,
                 generations=10,
                 nloci=7,
                 yaml_file='x'):

        if path.isfile(yaml_file):
            pass
            # TODO: populate all input parameters with input from yaml files
        else:
            self.popnumber = pop_number
            self.generations = generations
            self.nloci = nloci
            # TODO: check numvars is equal to length of minvalu and maxvalu
        self.adaptation_params = adaptation_params()

    def load_chemistry(self, chemistry="vesicle", chemistry_params=None, minvalu=None, maxvalu=None):

        builtin_chemistries = ["vesicle", "structure","micelle-structure","micelle-structure2",'structure-ml','structure-ml-train','structure-ml-train-one-comp','structure-binary','structure-ml-one','structure-ml-exp']
        self.chemistry = chemistry
        if chemistry in builtin_chemistries:
            sg = import_module('crease_ga_diameter.chemistries.' +
                               chemistry+'.scatterer_generator')
        else:
            raise CgaError('Currently unsupported shape {}'.format(chemistry))

        # TODO: Complete the checker
        if chemistry_params is None:
            self.scatterer_generator = sg.scatterer_generator()
        elif minvalu is None or maxvalu is None:
            warn("Unspecified minimum and/or maximum parameter boundaries. Fall back to the default minimum "
                 "and maximum parameter boundaries of shape {}.\n".format(chemistry), stacklevel=2)
            self.scatterer_generator = sg.scatterer_generator(chemistry_params)
            print("minimum parameter boundaries have been set to {},\n"
                  "maximum parameter boundaries have been set to {}.\n".format(
                      self.scatterer_generator.minvalu,
                      self.scatterer_generator.maxvalu))

        elif sg.scatterer_generator().numvars != len(minvalu) or sg.scatterer_generator().numvars != len(maxvalu):

            raise CgaError("Number of parameters in minvalu and/or maxvalu is not equal to number of parameters "
                           "required by shape {}.\n Shape {} requires {:d} parameters.\nminvalu has {:d} parameters.\n"
                           "maxvalu has {:d} parameters.".format(chemistry, chemistry, sg.scatterer_generator().numvars,
                                                                 len(minvalu), len(maxvalu)))
        else:
            self.scatterer_generator = sg.scatterer_generator(
                chemistry_params, minvalu, maxvalu)

        self.numvars = self.scatterer_generator.numvars
        self.minvalu = self.scatterer_generator.minvalu
        self.maxvalu = self.scatterer_generator.maxvalu

    def load_iq(self, input_file_path, q_bounds=None):
        loadvals = np.loadtxt(input_file_path)
        self.qrange_load = loadvals[:, 0]
        IQin_load = loadvals[:, 1]
        print('qrange start',self.qrange_load)
        if len(loadvals[0,:]) > 2:
            IQerr_load = loadvals[:,2]
        else:
            IQerr_load = None
        # self.IQin_load=np.true_divide(IQin_load,np.max(IQin_load))
        # TODO: highQ and lowQ needs to be able to be dynamically set
        if q_bounds is None:
            self.qrange = self.qrange_load
            self.IQin = IQin_load
            self.IQerr = IQerr_load
        else:
            lowQ = q_bounds[0]
            highQ = q_bounds[1]
            self.IQin = IQin_load[np.where(self.qrange_load >= lowQ)[
                0][0]:np.where(self.qrange_load <= highQ)[0][-1]+1]
            self.qrange = self.qrange_load[np.where(self.qrange_load >= lowQ)[
                0][0]:np.where(self.qrange_load <= highQ)[0][-1]+1]
            if IQerr_load is not None:
                self.IQerr = IQerr_load[np.where(self.qrange_load >= lowQ)[
                    0][0]:np.where(self.qrange_load <= highQ)[0][-1]+1]
            else:
                self.IQerr = IQerr_load

        baseline = self.IQin[0]
        self.IQin = np.true_divide(self.IQin, baseline)
        if IQerr_load is not None:
            self.IQerr = np.true_divide(self.IQerr,baseline)
        if self.chemistry == "structure":  
            self.scatterer_generator.setIQload(
                self.qrange, self.IQin, self.IQerr, self.popnumber, self.generations, self.nloci)
        if self.chemistry == "structure-binary":  
            IQin_load = loadvals[:, 1]
            IQin2_load = loadvals[:, 2]
            self.IQin = np.true_divide(IQin_load, IQin_load[0])
            self.IQin2 = np.true_divide(IQin2_load, IQin2_load[0])
            self.IQerr = None
            self.IQerr2 = None
            self.scatterer_generator.setIQload(
                self.qrange, self.IQin,self.IQin2, self.IQerr,self.IQerr2, self.popnumber, self.generations, self.nloci)
        if self.chemistry == "structure-ml-train":  
            self.scatterer_generator.setIQload(
                self.qrange, self.IQin, self.IQerr, self.popnumber, self.generations, self.nloci)
        if self.chemistry == "structure-ml-train-one-comp":  
            self.scatterer_generator.setIQload(
                self.qrange, self.IQin, self.IQerr, self.popnumber, self.generations, self.nloci)
        if self.chemistry == "structure-ml":  
            self.scatterer_generator.setIQload(
                self.qrange, self.IQin, self.IQerr, self.popnumber, self.generations, self.nloci)
        if self.chemistry == "structure-ml-one":  
            self.scatterer_generator.setIQload(
                self.qrange, self.IQin, self.IQerr, self.popnumber, self.generations, self.nloci)
        if self.chemistry == "structure-ml-exp":  
            self.scatterer_generator.setIQload(
                self.qrange, self.IQin, self.IQerr, self.popnumber, self.generations, self.nloci)
        if self.chemistry == "micelle-structure":  
            self.scatterer_generator.setIQload(
                self.qrange, self.IQin, self.IQerr, self.popnumber, self.generations, self.nloci)
        if self.chemistry == "micelle-structure2":  
            self.scatterer_generator.setIQload(
                self.qrange, self.IQin, self.IQerr, self.popnumber, self.generations, self.nloci)

    def solve(self, verbose=True, backend='debye', fitness_metric='chi2', output_dir='./', n_cores=1):
        ### checking if starting new run or restarting partial run
        if path.isfile('current_gen.txt'):
            currentgen = int(np.genfromtxt('current_gen.txt'))
            pop = np.genfromtxt('current_pop.txt')
            temp = np.genfromtxt('current_pm_pc.txt')
            pm = temp[0]
            pc = temp[1]
            self.adaptation_params.pc = pc
            self.adaptation_params.pm = pm
            print('Restarting from gen #{:d}'.format(currentgen+1))
        else:
            currentgen = 0
            pop = utils.initial_pop(self.popnumber, self.nloci, self.numvars)
            print('New run')
        self.cores = n_cores

        ## generate file to confirm supercomputer running correctly
        np.savetxt('supercomputer.txt',np.c_[0])

        for gen in range(currentgen, self.generations):
            startTime = time.time()
            if backend == 'debye':
                if self.chemistry == 'structure-ml-one' or self.chemistry == 'structure-ml-exp':
                    params = []
                    for val in range(self.popnumber):
                        # gets the current structure variables
                        param = utils.decode(
                            pop, val, self.nloci, self.minvalu, self.maxvalu)
                        params.append(param)
                    # calculate scattering
                    pool = mp.Pool(n_cores)
                    one = pool.starmap(self.scatterer_generator.doScatteringML, zip(
                        repeat(gen), range(self.popnumber), params))
                    pool.close()
                    pool.join()
                    one = np.array(one,dtype=object)
                    iq = one[:, 0]
                    iq = np.array(iq)
                    params = np.array(params)
                    pacc, gdm, elitei, IQid_str = self.scatterer_generator.fitness(pop, gen, output_dir, iq, params, metric=fitness_metric)
                    print('scattering time',(time.time()-startTime)/60.,'min')
                elif True:  
                    params = []
                    for val in range(self.popnumber):
                        # gets the current structure variables
                        param = utils.decode(
                            pop, val, self.nloci, self.minvalu, self.maxvalu)
                        params.append(param)
                        ## produce structures for LAMMPS
                        #self.scatterer_generator.produceStructure(param,val)
                    #"""
                    ## produce structures for LAMMPS
                    pool = mp.Pool(n_cores)
                    one = pool.starmap(self.scatterer_generator.produceStructure, 
                        zip(params,range(self.popnumber)))
                    pool.close()
                    pool.join()
                    print('datafiles created')
                    
                    # now run LAMMPS to create close-packed strucutres
                    run(['bash','run.sh',str(self.popnumber-1)])
                    while (int(check_output(['bash','run2.sh'])) < self.popnumber):
                        time.sleep(30)
                    #"""
                    #print('skip structure')
                    structureTime = time.time()
                    print('time for structure generation',(structureTime-startTime)/60.,'min')

                    # calculate scattering
                    pool = mp.Pool(n_cores)
                    one = pool.starmap(self.scatterer_generator.doScattering, zip(
                        repeat(gen), range(self.popnumber), params))
                    pool.close()
                    pool.join()
                    one = np.array(one,dtype=object)
                    atype = one[:, 0]
                    pos = one[:, 1]
                    iq = one[:, 2]
                    iq = np.array(iq)
                    if self.chemistry == "structure-ml-train":  
                        iq2 = one[:, 3]
                        iq2 = np.array(iq2)
                        print('saved iq values')
                        np.savetxt('genes_values_{}.txt'.format(gen),np.c_[params])
                        for qqq in range(self.popnumber):
                            with open('iqa_values_{}.txt'.format(gen),'a') as fiq:
                                np.savetxt(fiq,iq[qqq],newline=' ')
                                fiq.write('\n')
                            with open('iqb_values_{}.txt'.format(gen),'a') as fiq2:
                                np.savetxt(fiq2,iq2[qqq],newline=' ')
                                fiq2.write('\n')
                    if self.chemistry == "structure-ml-train-one-comp":  
                        fq = one[:,3]
                        sq = one[:,4]
                        print('saved iq values')
                        np.savetxt('genes_values_{}.txt'.format(gen),np.c_[params])
                        for qqq in range(self.popnumber):
                            with open('iq_values_{}.txt'.format(gen),'a') as fiq:
                                np.savetxt(fiq,iq[qqq],newline=' ')
                                fiq.write('\n')
                            with open('sq_values_{}.txt'.format(gen),'a') as fiq2:
                                np.savetxt(fiq2,sq[qqq],newline=' ')
                                fiq2.write('\n')
                            with open('fq_values_{}.txt'.format(gen),'a') as fiq3:
                                np.savetxt(fiq3,fq[qqq],newline=' ')
                                fiq3.write('\n')
                    
                    params = np.array(params)
                    if self.chemistry == 'structure-binary':
                        iq2 = np.array(one[:,3])
                        pacc, gdm, elitei, IQid_str = self.scatterer_generator.fitness(pop, gen, output_dir, atype, pos, iq,iq2, params, metric=fitness_metric)
                    else:
                        pacc, gdm, elitei, IQid_str = self.scatterer_generator.fitness(pop, gen, output_dir, atype, pos, iq, params, metric=fitness_metric)
                    print('scattering time',(time.time()-structureTime)/60.,'min')
                else:
                    pacc, gdm, elitei, IQid_str = self.fitness(
                        pop, gen, output_dir, metric=fitness_metric)
                    IQid_str = np.array(IQid_str)
            if self.chemistry == 'micelle-structure2':
                pop = self.scatterer_generator.genetic_operations(pop,gen,pacc,elitei,self.adaptation_params.pc,self.adaptation_params.pm)
            else:
                pop = self.genetic_operations(pop, pacc, elitei)
            self.adaptation_params.update(gdm)
            
            if self.chemistry == "structure-ml-train":  
                print('random pop for training ML')
                pop = utils.initial_pop(self.popnumber, self.nloci, self.numvars)
            if self.chemistry == "structure-ml-train-one-comp":  
                print('random pop for training ML')
                pop = utils.initial_pop(self.popnumber, self.nloci, self.numvars)
            
            ### save output from each generation
            np.savetxt('current_gen.txt',np.c_[gen])
            np.savetxt('current_pop.txt',np.c_[pop])
            np.savetxt('current_pm_pc.txt',np.c_[self.adaptation_params.pm,self.adaptation_params.pc])

            if verbose:
                figsize = (4, 4)
                fig, ax = plt.subplots(figsize=(figsize))
                ax.plot(self.qrange_load, self.IQin_load, color='k',
                        linestyle='-', ms=8, linewidth=1.3, marker='o')
                # ,marker='o')
                ax.plot(
                    self.qrange, IQid_str[elitei].transpose(), color='fuchsia', linestyle='-', ms=8, linewidth=2)
                plt.xlim(0.001, 0.1)
                plt.ylim(2*10**(-5), 20)
                plt.xlabel(r'q, $\AA^{-1}$', fontsize=20)
                plt.ylabel(r'$I$(q)', fontsize=20)
                ax.set_xscale("log")
                ax.set_yscale("log")
                plt.show()
        np.savetxt('simdone.txt',np.c_[gen])

    def fitness(self, pop, generation, output_dir, metric='log_sse'):
        tic = time.time()
        cs = 10
        F1 = open(output_dir+'z_temp_results_'+str(generation)+'.txt', 'w')
        np.savetxt(output_dir+'z_temp_population_' +
                   str(generation)+'.txt', np.c_[pop])

        fitn = np.zeros(self.popnumber)
        fitnfr = np.zeros(self.popnumber)
        fit = np.zeros(self.popnumber)
        qfin = self.qrange[-1]
        IQid_str = []
        params = []
        for val in range(self.popnumber):
            sys.stdout.write("\rGen {:d}/{:d}, individual {:d}/{:d}".format(
                generation+1, self.generations, val+1, self.popnumber))
            sys.stdout.flush()
            # gets the current structure variables
            param = utils.decode(pop, val, self.nloci,
                                 self.minvalu, self.maxvalu)
            params.append(param)

            ### calculate computed Icomp(q) ###
            IQid = self.scatterer_generator.converttoIQ(self.qrange, param)

            err = 0
            for qi, qval in enumerate(self.qrange):
                if (IQid[qi] > 0) & (self.IQin[qi] > 0):
                    if metric == 'log_sse':
                        if (qi < qfin):
                            # weighting factor
                            wil = np.log(np.true_divide(
                                self.qrange[qi+1], self.qrange[qi]))
                        else:
                            # weighting factor
                            wil = np.log(np.true_divide(
                                self.qrange[qi], self.qrange[qi-1]))
                        # squared log error
                        err += wil * \
                            (np.log(np.true_divide(
                                self.IQin[qi], IQid[qi])))**2
                    elif metric == 'ch2':
                        # chi^2 with weighting of sqrt(IQin)
                        err += np.true_divide(
                            np.square(self.IQin[qi]-IQid[qi]), self.IQin[qi])

            fit[val] = err
            IQid_str.append(IQid)

            F1.write(
                (str(val)+' '+str(param[0])+' '+str(param[1])+' '+str(param[2])+str(err)+'\n'))
            F1.flush()

        maxerr = np.max(fit)  # determines maximum SSerror for the population
        fitn = fit-maxerr  # determines error differences
        # maxfit=np.max(fitn)
        sumup = np.sum(fitn)

        avgfit = np.true_divide(sumup, self.popnumber)
        dval = -avgfit
        # linear scaling with cs as a scaleFactor
        ascale = np.true_divide(avgfit, dval)*(cs-1.0)
        bscale = avgfit*(1.0-ascale)

        sumup = 0

        # get scaled fitness to enable selection of bad candidates
        for val in range(self.popnumber):
            if (fitn[val] > avgfit):
                fitnfr[val] = ascale*fitn[val]+bscale
            else:
                fitnfr[val] = fitn[val]

        sumup = np.sum(fitnfr)

        pacc = np.zeros(self.popnumber)
        prob = np.true_divide(fitnfr, sumup)
        pacc = np.cumsum(prob)

        ### returns cummulative relative error from which individuals can be selected ###
        maxfit = np.min(fit)
        elitei = np.where(fit == maxfit)[0]                  # Best candidate
        secondfit = sorted(fit)[1]
        # Second best candidate
        secondi = np.where(fit == secondfit)[0]
        avgfit = np.average(fit)
        avgi = np.array([(np.abs(fit-avgfit)).argmin()])   # Average candidate
        minfit = np.max(fit)
        mini = np.where(fit == minfit)[0]                    # Worst candidate
        if avgfit == 0:
            avgfit = 1
        gdm = np.true_divide(maxfit, avgfit)
        if len(elitei) > 1:
            elitei = elitei[0]
        if len(secondi) > 1:
            secondi = secondi[0]
        if len(avgi) > 1:
            avgi = avgi[0]
        if len(mini) > 1:
            mini = mini[0]

        f = open(output_dir+'fitness_vs_gen.txt', 'a')
        if generation == 0:
            f.write('gen mini min avgi avg secondi second besti best\n')
        f.write('%d ' % (generation))
        f.write('%d %.8lf ' % (mini, minfit))
        f.write('%d %.8lf ' % (avgi, avgfit))
        f.write('%d %.8lf ' % (secondi, secondfit))
        f.write('%d %.8lf ' % (elitei, maxfit))
        f.write('\n')
        f.close()
        print('Generation time: {:.3f}s'.format(time.time()-tic))
        print('Generation best fitness: {:.4f}'.format(maxfit))
        print('Generation best fitness: {:.3f}'.format(gdm))
        params = np.array(params)
        print('Generation best parameters '+str(params[elitei]))

        return pacc, gdm, elitei, IQid_str

    def genetic_operations(self, pop, pacc, elitei):
        popn = np.zeros(np.shape(pop))
        cross = 0
        mute = 0
        pc = self.adaptation_params.pc
        pm = self.adaptation_params.pm

        for i in range(self.popnumber-1):

            #####################    Crossover    ####################

            # Selection based on fitness
            testoff = random.random()
            isit = 0
            npart1 = 1
            for j in range(1, self.popnumber):
                if (testoff > pacc[j-1]) & (testoff < pacc[j]):
                    npart1 = j

            testoff = random.random()
            isit = 0
            npart2 = 1
            for j in range(self.popnumber):
                if (testoff >= pacc[j-1]) & (testoff != pacc[j]):
                    npart2 = j

            # Fit parents put in array popn

            popn[i, :] = pop[npart1, :]

            testoff = random.random()
            loc = int((testoff*(self.numvars-1))*self.nloci)
            if loc == 0:
                loc = self.nloci
            testoff = random.random()

            # crossover
            if (testoff <= pc):
                cross += 1
                popn[i, loc:] = pop[npart2, loc:]

        #####################    Mutation    ####################

            for j in range(self.nloci*self.numvars):
                testoff = random.random()
                if (testoff <= pm):
                    popn[i, j] = random.randint(0, 1)
                    mute += 1

        #####################    Elitism    ####################

        popn[-1, :] = pop[elitei, :]

        print('pc', pc)
        print('#crossovers', cross)
        print('pm', pm)
        print('#mutations', mute)
        print('\n')

        return popn
