#calculates omega or intra-micellar structure factor
import numpy as np #allows efficient array calculations
import numexpr as ne #This allows for speedup of the calculation

def LPOmega(qrange, nAin, nAout, nB, r):                # qvalues number_of_B number_of_A scatterer_coordinates
    Ntot=nAin+nB+nAout                                  # Total number of scatterers to loop through
    omegaarrt=np.zeros((1,len(qrange)))                 # initiating array
    for step in range(0, 1):                            # loops through independent_scatterer_placemetns
        omegaarr=np.zeros((1,len(qrange)))              # initiating array
        rur=r[step,:,:]                                 # selects      
        for i in range(Ntot-1):                         # loops through index and all further indexes to prevent double counting 
            x = np.square(rur[0,i]-rur[0,(i+1):])
            y = np.square(rur[1,i]-rur[1,(i+1):])
            z = np.square(rur[2,i]-rur[2,(i+1):])
            rij = np.sqrt(np.sum([x,y,z],axis=0))       # calculates the distances
            rs = rij[:,np.newaxis]                      # reshapes array for consistency
            Q = qrange[np.newaxis,:]                    # reshapes array for consistency
            vals = ne.evaluate("sin(Q*rs)/(Q*rs)")      # ne is efficient at calculations
            inds=np.argwhere(np.isnan(vals))            # error catching in case there are NaN values
            if len(inds)>0:
                for val in inds:
                    vals[val[0],val[1]]=1
            inds_double_check=np.argwhere(np.isnan(vals))
            if len(inds_double_check)>0:
                print('nan error!')
            vals = ne.evaluate("sum((vals), axis=0)")   # adds together scatterer contributions for each q value
            omegaarr+=vals
        omegaarr=np.true_divide(2*omegaarr,Ntot)+1      # 1 accounts for the guarenteed overlap of same bead  # 2* accounts for double counting avoided to reduce computational expense by looping for all other pairs
        omegaarrt+=omegaarr                             # stores values between loops
#     omegaarrt=np.true_divide(omegaarrt,nLP)             # averages between independent_scatterer_placements
#     omegaarrt=omegaarrt.reshape(len(qrange),)           # reshapes array for output
    return omegaarrt