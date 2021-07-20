import numpy as np
from Initial_pop import Initial_pop
from decode import decode
from fitness import fitness
import random    
import matplotlib
matplotlib.use('Agg') ## uncomment this when running on cluster, comment out this line if on local
import matplotlib.pyplot as plt
import sys

##### Initial Parameters #####
## comments
gdmmin=0.005
gdmmax=0.85
pcmin=0.5
pcmax=1
pmmin=0.006
pmmax=0.25
kgdm=1.1
pc=0.6
pm=0.001

##################### generate initial population############################## 
popnumber=80 #population number
generations=100 # number of generations
#popnumber=int(sys.argv[1]) #population number
#generations=int(sys.argv[2]) # number of generations
print(popnumber, generations)

nloci=7 #number of binary numbers per variable, 7 for each variable (ranged from 0 to 127) 

#### Min max values of parameters ####
minvalu=(50, 30, 30, 30, 0.1, 0.0, 0)    # Rcore dR_Ain dR_B dR_Aout split_A sigmaR Background
maxvalu=(400, 200, 200, 200, 0.45, 0.45, 4)    # Rcore dR_Ain dR_B dR_Aout split_A sigmaR Background

minvalu=np.array(minvalu)
maxvalu=np.array(maxvalu)
numvars=len(minvalu) #number of variables 
pop=Initial_pop(popnumber, numvars, nloci) # #inviduals #variables #parameters


##################### set up q and read in IQ ############################## 
loadvals=np.loadtxt('Itot_disper_10_Ain12_B6_Aout12_nLP7_dR0.2.txt')
qrange_load=loadvals[:,0]
IQin_load=loadvals[:,1]
IQin_load=np.true_divide(IQin_load,np.max(IQin_load))

num_scatterers=54  # in JACS paper, using 8 scatterers

N = 54                   # Number of beads on chain
rho_B = 0.5              # density/volume fraction of LJ bead in B layer
MB = np.pi/6*(50.4)**3   # volume of a B monomer
lmono_a = 50.4           # Angstrom 'monomer contour length' for simulations 50.4
lmono_b = 50.4           # Angstrom 'monomer contour length' for simulations 50.4
sigmabead = np.true_divide(N*lmono_b,num_scatterers)
fb = 0.6                 ## Keep constant at 0.6
nLP = 7

set_lowq=0.0049 ## 0.0049 corresponds to 1282 Ang
set_highq=0.07 ##0.078

highQ=qrange_load[np.where(qrange_load<=set_highq)[0]][-1]
lowQ=qrange_load[np.where(qrange_load>=set_lowq)[0]-1][0]
print(highQ)
print(lowQ)

Results=[]
Results.append(popnumber)
Results.append(nloci)
Results.append(num_scatterers)
Results.append(N)
Results.append(fb)
Results.append(rho_B)
Results.append(MB)
Results.append(lmono_a)
Results.append(lmono_b)
Results.append(sigmabead)
Results.append(set_lowq)
Results.append(set_highq)
Results.append(highQ)
Results.append(lowQ)
Results=np.array(Results)
np.savetxt('details.txt',np.c_[Results])

Results=[]
Results.append(minvalu)
Results.append(maxvalu)
np.savetxt('val_bounds.txt',np.c_[Results])

qrange = qrange_load[ np.where(qrange_load==lowQ)[0][0]:np.where(qrange_load==highQ)[0][0] ]
IQin = IQin_load[ np.where(qrange_load==lowQ)[0][0]:np.where(qrange_load==highQ)[0][0] ]

baseline = IQin[0]
IQin = np.true_divide(IQin,baseline)
IQin_load = np.true_divide(IQin_load,baseline)
generation = 0
pacc, gdm, elitei, secondi, maxfit, secondfit, avgfit, minfit, IQid_str, avgi, mini, params = fitness(
    IQin, qrange, pop, nloci, N, fb, sigmabead, minvalu, maxvalu, rho_B, MB, lmono_a, lmono_b, 
    num_scatterers, nLP, generation
    )
 
f = open( 'fitness_vs_gen.txt', 'a' )
if generation == 0:
    f.write( 'gen min avg second best\n' )
f.write( '%d ' %(generation) )
f.write( '%.8lf ' %minfit )
f.write( '%.8lf ' %avgfit )
f.write( '%.8lf ' %secondfit )
f.write( '%.8lf ' %maxfit )
f.write( '\n' )
f.close()
str_params = []
str_params.append(params)
str_IQid_str = []
str_IQid_str.append(IQid_str)

for gen in range(generations):

    popn=np.zeros(np.shape(pop))

    if (gdm>gdmmax):
        pm=pm*kgdm
        pc=np.true_divide(pc,kgdm)
    elif (gdm<gdmmin):
        pm=np.true_divide(pm,kgdm)
        pc=pc*kgdm
    if (pm>pmmax):
        pm=pmmax
    if (pm<pmmin):
        pm=pmmin
    if (pc>pcmax):
        pc=pcmax
    if (pc<pcmin):
        pc=pcmin
    cross=0
    mute=0
    
    
    for i in range(popnumber-1):
        
        #####################    Crossover    ####################
        
        #Selection based on fitness
        testoff=random.random()
        isit=0
        npart1=1
        for j in range(1,popnumber):
            if (testoff>pacc[j-1])&(testoff<pacc[j]):
                npart1=j

        testoff=random.random()
        isit=0
        npart2=1
        for j in range(popnumber):
            if (testoff>=pacc[j-1])&(testoff!=pacc[j]):
                npart2=j

        #Fit parents put in array popn
        
        popn[i,:]=pop[npart1,:]
 
                
        testoff=random.random()
        loc=int((testoff*(numvars-1))*nloci)
        if loc==0:
            loc=nloci
        testoff=random.random()

        #crossover
        if (testoff<=pc):
            cross+=1
            popn[i,loc:]=pop[npart2,loc:]


    #####################    Mutation    ####################


        for j in range(nloci*numvars):
            testoff=random.random()
            if (testoff<=pm):
                popn[i,j]=random.randint(0,1)
                mute+=1


    #####################    Elitism    ####################
    
    popn[-1,:]=pop[elitei,:]

    pop=popn
    print('gen',gen+1)
    print('pc',pc)
    print('#crossovers',cross)
    print('pm',pm)
    print('#mutations',mute)
    print('gdm',gdm) 
    pacc, gdm, elitei, secondi, maxfit, secondfit, avgfit, minfit, IQid_str, avgi, mini, params = fitness(
        IQin, qrange, pop, nloci, N, fb, sigmabead, minvalu, maxvalu, rho_B, MB, lmono_a, lmono_b, 
        num_scatterers, nLP, generation
        )
    print('maxfit2',maxfit)
    print('                   ')

    np.savetxt('zzz_pop_gen'+str(gen)+'.txt',np.c_[pop])
    val_str=[]
    val_str.append(pc)
    val_str.append(pm)
    val_str.append(gdm)
    val_str.append(elitei)
    val_str=np.array(val_str)
    np.savetxt('zzz_pc_pm_gdm.txt',np.c_[val_str])
    
    np.savetxt('current_pop.txt',np.c_[pop])
    np.savetxt('current_pc_pm_gdm.txt',np.c_[val_str])
    np.savetxt('current_gen.txt',np.c_[gen])
    np.savetxt('current_pacc.txt',np.c_[pacc])
    
    f = open( 'fitness_vs_gen.txt', 'a' )
    f.write( '%d ' %(gen+1) )
    f.write( '%.8lf ' %minfit )
    f.write( '%.8lf ' %avgfit )
    f.write( '%.8lf ' %secondfit )
    f.write( '%.8lf ' %maxfit )
    f.write( '\n' )
    f.close()
    
    np.savetxt('done.txt',np.c_[(1)])
    
    np.savetxt('parameters0.txt',np.c_[params[0]])#,fmt='%i') 
    np.savetxt('parametersfin.txt',np.c_[params[-1]])#,fmt='%i') 

    elitei=int(elitei)
    secondi=int(secondi)
    print(avgfit, avgi, mini)
    avgi=int(avgi)
    print('test1')
    mini=int(mini)
    print('test2') 
    
    print('opt params:')
    print('Rcore', np.round(params[elitei][0]/10,1), 'nm')
    print('dR_Ain', np.round(params[elitei][1]/10,1), 'nm')
    print('dR_B', np.round(params[elitei][2]/10,1), 'nm')
    print('dR_Aout', np.round(params[elitei][3]/10,1), 'nm')
    print('splitA', np.round(params[elitei][4],1))
    print('sigmaR', np.round(params[elitei][5],1))
    print('Background exponent',np.round(-1*params[elitei][5],1))
    np.savetxt('opt_Itot.txt',np.c_[IQid_str[elitei]])
    
    Results=[]
    Results.append(params[elitei][0]/10)
    Results.append(params[elitei][1]/10)
    Results.append(params[elitei][2]/10)
    Results.append(params[elitei][3]/10)
    Results.append(params[elitei][4])
    Results.append(params[elitei][5])
    Results.append(np.round(-1*params[elitei][5],1))
    np.savetxt('Results_fin.txt',np.c_[Results])
    
    ### save optimal parameters each generation ###
    f = open( 'opt_parameters.txt', 'a' )
    f.write( '%d ' %gen )
    for i in range( 0, len(params[0]) ):
        f.write( '%.6lf ' %params[elitei][i] )
    f.write( '\n' )
    f.close()
    
    f = open( 'second_parameters.txt', 'a' )
    f.write( '%d ' %gen )
    for i in range( 0, len(params[0]) ):
        f.write( '%.6lf ' %params[secondi][i] )
    f.write( '\n' )
    f.close()

    f = open( 'avg_parameters.txt', 'a' )
    f.write( '%d ' %gen )
    for i in range( 0, len(params[0]) ):
        f.write( '%.6lf ' %params[avgi][i] )
    f.write( '\n' )
    f.close()

    f = open( 'min_parameters.txt', 'a' )
    f.write( '%d ' %gen )
    for i in range( 0, len(params[0]) ):
        f.write( '%.6lf ' %params[mini][i] )
    f.write( '\n' )
    f.close()

    np.savetxt('IQelite_gen%d.txt' %gen, np.c_[qrange, IQid_str[elitei]], fmt='%0.8f')
    np.savetxt('IQsecond_gen%d.txt' %gen, np.c_[qrange, IQid_str[secondi]], fmt='%0.8f')
    np.savetxt('IQavg_gen%d.txt' %gen, np.c_[qrange, IQid_str[avgi]], fmt='%0.8f')
    np.savetxt('IQmin_gen%d.txt' %gen, np.c_[qrange, IQid_str[mini]], fmt='%0.8f')
    
    figsize=(4,4)
    fig, ax = plt.subplots(figsize=(figsize))
    ax.plot(qrange_load,IQin_load,color='k',linestyle='-',ms=8,linewidth=1.3,marker='o')
    ax.plot(qrange,IQid_str[elitei],color='fuchsia',linestyle='-',ms=8,linewidth=2)#,marker='o')
    plt.xlim(0.001,0.1)
    plt.ylim(2*10**(-5),20)
    plt.xlabel(r'q, $\AA^{-1}$',fontsize=20)
    plt.ylabel(r'$I$(q)',fontsize=20)
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.show()

