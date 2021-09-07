import numpy as np
def initial_pop(popnumber, nloci, numvars):
    '''
    Produce a generation of (binary) chromosomes.
    
    Parameters
    ----------
    popnumber: int
        Number of individuals in a population.
    nloci: int
        Number of binary bits to represent each parameter in a chromosome.
    numvars: int
        Number of parameters in a chromosome.
        
    Returns
    -------
    pop: np.array of size (`popnumber`,`nloci`*`numvars`)
        A numpy array of binary bits representing the entire generation, 
        with each row representing a chromosome.
    '''
    pop=np.zeros((popnumber,nloci*numvars))
    for i in range(popnumber):
        for j in range(nloci*numvars):
            randbinary=np.random.randint(2)
            pop[i,j]=randbinary
    return pop
