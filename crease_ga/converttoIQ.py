# calculates the I(q) for the individual in the population
from decode import decode
from genLP import genLP
from LPFbead import LPFbead
from LPOmega import LPOmega
import numpy as np

def converttoIQ(qrange, param, sigmabead, N, fb, rho_B, MB, lmono_a, lmono_b, 
                num_scatterers, nLP): 
    ##q values, decoded parameters, 
    # number of repeat units per chain, fraction of B beads per chain, core density, 
    # scatterer diameter, molar mass of B chemistry, 
    # length of A chemistry bond, length of B chemistry bond, 
    # number of scatterers per chain, # of replicates, stdev in Rcore size
    
    IQid=np.zeros((len(qrange)))      #initiates array for output IQ
    
    #print(param)
    ######### Calculate inputs to generating scatterer placements ########
    Rcore=param[0]
    dR_Ain=param[1]
    dR_B=param[2]
    dR_Aout=param[3]
    sAin=param[4]     # split of type A scatterer 
    sigmaR=param[5]   # variation in Rcore, dispersity
    #print(Rcore, dR_Ain, dR_B, dR_Aout, sAin)
    Background=10**(-param[6]) 
    
    varR = Rcore*sigmaR # variation in Rcore
    disper = np.array([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0])
    sum_omegaarr=np.zeros((1,len(qrange)))
    
    for step in range(0, nLP):
        Rcore = param[0] + varR*disper[step + int((9-nLP)/2.)]          ## add displacement to Rcore
        print("disper = ", disper[step + int((9-nLP)/2.)])
        print("Rcore = ", Rcore)
        vol_B = (4/3.0)*np.pi*(np.power(Rcore + dR_Ain + dR_B, 3) 
                               - np.power(Rcore + dR_Ain, 3)) ## volume of solvophobic layer B
        nagg = int(np.true_divide( rho_B*vol_B, N*fb*MB ))    ## number of chains in vesicle
        ntot = nagg*num_scatterers                            ## total number of scatterers
        nB = int(ntot*fb)                                     ## number of scatterers in B
        nAin = int(ntot*(1-fb)*sAin)                          ## number of scatterers in A_in
        nAout = int(ntot*(1-fb)*(1-sAin))                     ## number of scatterers in A_out
        #print(vol_B, nagg, ntot, nB, nAin, nAout)

        for reps in range(0, 1):
            ### Generates scatterer positions in structure ###
            r = genLP(Rcore, dR_Ain, dR_B, dR_Aout, sigmabead, nAin, nAout, nB)

            ### Calculates omega from scatterers in shape ###
            sum_omegaarr += LPOmega(qrange, nAin, nAout, nB, r)
    
    omegaarr=np.true_divide(sum_omegaarr,nLP*1)          # average omega
    omegaarr=omegaarr.reshape(len(qrange),)
    Fqb=LPFbead(qrange,sigmabead)                      # calcualtes sphere shape factor
    F2qb=np.multiply(Fqb,Fqb)                          # Sphere shape factor square
    sqmm=np.ones((np.shape(Fqb)))                      # assuming dilute mixture the micelle-micelle structure factor = 1
    F2qb_sqmm=np.multiply(F2qb,sqmm)                   # determines the micelle form factor
    IQid=np.multiply(omegaarr,F2qb_sqmm)+Background    # calculates Icomp
    maxIQ=np.max(IQid)                                  
    IQid=np.true_divide(IQid,maxIQ)                    # normalizes the I(q) to have its maximum = 1
    return IQid
