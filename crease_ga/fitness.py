#this code determines the fitness of all the individuals in the GA
from decode import decode
from converttoIQ import converttoIQ
import numpy as np
import time

'''
notes
popnumber = number of individuals
pop stores pop and genes
nloci = 7* #variables for the binary representation
'''
def fitness(IQin, qrange, pop, nloci, N, fb, sigmabead, minvalu, maxvalu, 
            rho_B, MB, lmono_a, lmono_b, num_scatterers, nLP, generation):
    cs=10
    F1= open('z_temp_results_'+str(generation)+'.txt','w')
    np.savetxt('z_temp_population_'+str(generation)+'.txt',np.c_[pop])
    popnumber=len(pop[:,0])
    fitn=np.zeros(popnumber)
    fitnfr=np.zeros(popnumber)
    fit=np.zeros(popnumber)
    qfin=qrange[-1]
    IQid_str=[]
    params=[]
    for val in range(popnumber):
        t0=time.time()
        param=decode(pop, val, nloci, minvalu, maxvalu) #gets the current micelle variables
        params.append(param)

        IQid=converttoIQ(qrange, param, sigmabead, N, fb, rho_B, MB, lmono_a, lmono_b, 
                         num_scatterers, nLP)
        
        err=0
        for qi,qval in enumerate(qrange):
            if (IQid[qi]>0)&(IQin[qi]>0):
                if (qi<qfin):
                    wil=np.log(np.true_divide(qrange[qi+1],qrange[qi]))  #weighting factor
                else:
                    wil=np.log(np.true_divide(qrange[qi],qrange[qi-1]))  #weighting factor
                err+=wil*(np.log(np.true_divide(IQin[qi],IQid[qi])))**2  #squared log error 
        fit[val]=err
        IQid_str.append(IQid)
        
        t1=time.time()
        dt=t1-t0
        
        F1.write((str(val)+' '+str(dt)+' '+str(param[0])+' '+str(param[1])+' '+str(param[2])+str(err)+'\n'))
        F1.flush()

            
    maxerr=np.max(fit)      #determines maximum SSerror for the population
    fitn=np.subtract(maxerr,fit)#determines error differences?
    maxfit=np.max(fitn)
    sumup=np.sum(fitn)
    
    avgfit=np.true_divide(sumup,popnumber)
    dval=maxfit-avgfit
    ascale=np.true_divide(avgfit,dval)*(cs-1.0)     #linear scaling with cs as a scaleFactor
    #max - avg*sc
    bscale=avgfit*(1.0-ascale)

    sumup=0
    
    #get scaled fitness to enable selection of bad candidates
    for val in range(popnumber):
        if (fitn[val]>avgfit):
            fitnfr[val]=ascale*fitn[val]+bscale
        else:
            fitnfr[val]=fitn[val]

    sumup=np.sum(fitnfr)

    pacc=np.zeros(popnumber)
    prob=np.true_divide(fitnfr,sumup)
    pacc=np.cumsum(prob)

    # returns cummulative relative error from which individuals can be selected
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
    return pacc, gdm, elitei, secondi, maxfit, secondfit, avgfit, minfit, IQid_str, avgi, mini, params