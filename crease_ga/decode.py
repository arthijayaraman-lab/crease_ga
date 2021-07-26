# changes binary representation of variables to values between set min and max of each variable
import numpy as np

def decode(pop, indiv, nloci, minvalu, maxvalu):
    import numpy as np
    
    deltavalu=maxvalu-minvalu
    nvars=len(minvalu)
    valdec=np.zeros(nvars)
    param=np.zeros(nvars)
    #   decodes from binary to values between max and min
    for j in range(nvars): 
        n=nloci
        for i in range(j*nloci,(j+1)*nloci):
            n=n-1
            valdec[j]+=pop[indiv,i]*(2**n)
            
        param[j]=minvalu[j]+np.true_divide((deltavalu[j])*(valdec[j]),2**nloci)
    return param