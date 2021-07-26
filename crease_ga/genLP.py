#this generates placement of scatterers for spherical micelles
import numpy as np
import random   

def gen_layer(rin, rout, nsize):
    R = 1.0
    
    phi = np.random.uniform(0, 2*np.pi, size=(nsize))
    costheta = np.random.uniform(-1, 1, size=(nsize))
    u = np.random.uniform(rin**3, rout**3, size=(nsize))

    theta = np.arccos( costheta )
    r = R * np.cbrt( u )
    
    x = r * np.sin( theta ) * np.cos( phi )
    y = r * np.sin( theta ) * np.sin( phi )
    z = r * np.cos( theta )
    
    return( x, y, z )

def genLP(Rcore, dR_Ain, dR_B, dR_Aout, sigmabead, nAin, nAout, nB):  
    # core radius, inner A layer thickness, B layer thickness, outer A layer thickness, 
    # bead diameter, # of inner A beads, # of outer A beads, # of B beads
    
    ntot = nAin+nB+nAout
    print('ntot is %d' %ntot)
    power = 2
    r = np.zeros((1, 3, ntot))
    types = np.zeros((ntot))
    
    ### Create configuration for each replicate with dispersity ###
    for step in range(0, 1):
        ### Populate A inner Layer ###
        x, y, z = gen_layer(Rcore, Rcore+dR_Ain, nAin)
        for i in range(nAin):
            r[0,:,i] = [x[i], y[i], z[i]]
            types[i] = 1

        ### Populate B middle Layer ###
        x, y, z = gen_layer(Rcore+dR_Ain, Rcore+dR_Ain+dR_B, nB)
        for i in range(nB):
            r[0,:,i+nAin] = [x[i], y[i], z[i]]
            types[i+nAin] = 2

        ### Populate A outer Layer ###
        x, y, z = gen_layer(Rcore+dR_Ain+dR_B, Rcore+dR_Ain+dR_B+dR_Aout, nAout)
        for i in range(nAout):
            r[0,:,i+nAin+nB] = [x[i], y[i], z[i]]
            types[i+nAin+nB] = 1            
    return r
