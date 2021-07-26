import numpy as np
def initial_pop(popnumber, nloci, numvars):
        
    pop=np.zeros((popnumber,nloci*numvars))
    for i in range(popnumber):
        for j in range(nloci*numvars):
            randbinary=random.randint(0,1)
            pop[i,j]=randbinary
    return pop