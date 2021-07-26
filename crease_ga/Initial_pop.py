#This determines the inital random binary sequences for the geness of the whole population
'''
notes
popnumber = number of individuals
pop stores pop and genes
if nloci = 7 then 7 * #variables for the binary representation
'''
import random    
import numpy as np
    
def Initial_pop(popnumber, numvars, nloci):
    pop=np.zeros((popnumber,nloci*numvars))
    for i in range(popnumber):
        for j in range(nloci*numvars):
            randbinary=random.randint(0,1)
            pop[i,j]=randbinary
    return pop