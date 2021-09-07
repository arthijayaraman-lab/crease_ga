import numpy as np
def decode(pop, indiv, nloci, minvalu, maxvalu):
    '''
    Convert a binary chromosome from a generation back to decimal parameter values.

    Parameters
    ----------
    pop: np.array.
        A numpy array of binary bits representing the entire generation, 
        with each row representing a chromosome.
    indiv: int.
        The row ID of the chromosome of interest.
    nloci: int
        Number of binary bits used to represent each parameter in a chromosome.
    minvalu, maxvalu: list-like.
        The minimum/maximum boundaries (in decimal value) of each parameter.
        "All-0s" in binary form will be converted to the minimum for a
        parameter, and "all-1s" will be converted to the maximum.

    Returns
    -------
    param: np.array.
        A 1D array containing the decimal values of the input parameters.
    '''
    import numpy as np
    minvalu = np.array(minvalu)
    maxvalu = np.array(maxvalu)
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
