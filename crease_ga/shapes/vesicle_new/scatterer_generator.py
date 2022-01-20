import numpy as np
import random
import numexpr as ne

def gen_layer(rin, rout, nsize):
        R = 1.0

        phi = np.random.uniform(0, 2*np.pi, size=(nsize))
        costheta = np.random.uniform(-1, 1, size=(nsize))
        u = np.random.uniform(rin**3, rout**3, size=(nsize)) 

        theta = np.arccos( costheta )
        r = R * np.cbrt( u ) #ensures an even distribution

        x = r * np.sin( theta ) * np.cos( phi )
        y = r * np.sin( theta ) * np.sin( phi )
        z = r * np.cos( theta )
        
        ret = np.zeros((len(x),3))
        ret[:,0] = x
        ret[:,1] = y
        ret[:,2] = z
        return ret
    
def LPFbead(qrange, sigmabead):
    '''
    Compute the spherical form factor given a range of q values.
    
    Parameters
    ----------
    qrange: numpy.array
        array of values in q-space to compute form factor for.
    sigmabead: float
        diameter of the sphere.
    
    Return
    -------
    Fqb: numpy.array
        array of values of the spherical form factors (F(q)) computed at q-points listed in qrange.
    '''
    
    R=np.true_divide(sigmabead,2)
    QR=np.multiply(qrange,R)
    Fqb=np.multiply(np.true_divide(np.sin(QR)-np.multiply(QR,np.cos(QR)),np.power(QR,3)),3)  

    return Fqb

def LPOmega(qrange, r, sigmabead, sld_Ain, sld_B, sld_Aout):                # qvalues number_of_B number_of_A scatterer_coordinates
    
    omegaarr=np.zeros(len(qrange))              # initiating array
    sff = LPFbead(qrange,sigmabead)
    sld = [sld_Ain, sld_B, sld_Aout]
    for ri in range(len(r)):
        rur = r[ri]
        fi = sff*sld[ri]
        for i in range(len(rur)):
            all_disp = rur[i,:]-rur[(i+1):,:]
            rij = np.sqrt(np.sum(np.square(all_disp),axis=1))
            rij = rij.transpose()
            rs = rij[:,np.newaxis]                      # reshapes array for consistency
            Q = qrange[np.newaxis,:]                    # reshapes array for consistency
            vals = ne.evaluate("sin(Q*rs)/(Q*rs)")      # ne is efficient at calculations
            inds=np.argwhere(np.isnan(vals))            # error catching in case there are NaN values
            if len(inds)>0:
                for val in inds:
                    vals[val[0],val[1]]=1
            vals = ne.evaluate("sum((vals), axis=0)")   # adds together scatterer contributions for each q value
            omegaarr += fi**2*(2*vals+1)      # 1 accounts for the guarenteed overlap of same bead  # 2* accounts for double counting avoided to reduce computational expense by looping for all other pairs
            for rj in range(ri+1,len(r)):
                rsec = r[rj]
                all_disp = rur[i,:]-rsec
                rij = np.sqrt(np.sum(np.square(all_disp),axis=1))
                rij = rij.transpose()
                rs = rij[:,np.newaxis]                      # reshapes array for consistency
                Q = qrange[np.newaxis,:]                    # reshapes array for consistency
                vals = ne.evaluate("sin(Q*rs)/(Q*rs)")      # ne is efficient at calculations
                inds=np.argwhere(np.isnan(vals))            # error catching in case there are NaN values
                if len(inds)>0:
                    for val in inds:
                        vals[val[0],val[1]]=1
                vals = ne.evaluate("sum((vals), axis=0)")   # adds together scatterer contributions for each q value
                fj = sff*sld[rj]
                omegaarr += fi*fj*vals
            

    return omegaarr

def visualize(r, Rcore, dR_Ain, dR_B, dR_Aout, sigmabead):
    import py3Dmol
    view = py3Dmol.view()
    
    for rur in r:
        for ri in rur:
            if np.linalg.norm(ri) < Rcore+dR_Ain or np.linalg.norm(ri) > (Rcore+dR_Ain+dR_B):
                col = 'blue'
            else:
                col = 'red'
            view.addSphere(
                {
                    'center': {'x': ri[0], 'y': ri[1], 'z': ri[2]},
                           'radius': sigmabead/2,
                           'color': col,
                           'alpha': 0.9,
                }
                          )
        #view.zoomTo()
        view.show()

        return view
    
           

def genLP(R_core, dR_Ain, dR_B, dR_Aout, scat_density):  
        # core radius, inner A layer thickness, B layer thickness, outer A layer thickness, 
        # bead diameter, # of inner A beads, # of outer A beads, # of B beads

        r = [[],[],[]]
        ### Create configuration for each replicate with dispersity ###
       
        nAin = int(4/3*np.pi*((R_core+dR_Ain)**3-R_core**3)*scat_density)
        r[0] = gen_layer(R_core,R_core+dR_Ain,nAin)
        nB = int(4/3*np.pi*((R_core+dR_Ain+dR_B)**3-(R_core+dR_Ain)**3)*scat_density)
        r[1] = gen_layer(R_core+dR_Ain,R_core+dR_Ain+dR_B,nB)
        nAout = int(4/3*np.pi*((R_core+dR_Ain+dR_B+dR_Aout)**3-(R_core+dR_Ain+dR_B)**3)*scat_density)
        r[2] = gen_layer(R_core+dR_Ain+dR_B,R_core+dR_Ain+dR_B+dR_Aout,nAout)
        
        return r

    
class scatterer_generator:
    '''
    The wrapper class for vesicle shape. Default length unit: Angstrom.
    

    Notes
    -----
    **The following 7 shape-specific descriptors are to be specified by user (see
    *Attributes*) as 
    a list, in the precise order as listed, while calling `Model.load_shape`
    to load this shape:**


    eta_B:
        scatterer density. Default: 0.5
    sigmabead:
        scatterer diameter. Default: 1


    **The following 7 parameters are to be predicted, in the precise order
    as listed, by GA:**
    
    R_total:
        R_core+R_in+R_b+R_out. Default [min,max]: [20 A, 5000 A]
    f_core:
        R_core/R_total. Default [min,max]: [0.1,0.95]
    f_Ain:
        R_Ain/(R_total-R_core). Default [min,max]: [0.01,0.99]
    f_Aout:
        R_Aout/(R_total-R_core-R_Ain). Default [min,max]: [0.01,0.99]
    sld_Ain:
        SLD of Ain layer. 
        Default [min,max]: [-1,1]
    sld_B:
        SLD of B layer. 
        Default [min,max]: [-1,1]
    sld_Aout:
        SLD of Aout layer. 
        Default [min,max]: [-1,1]

    log10(bg):
        Negative log10 of background intensity. 
        E.g. an background intensity of 0.001 leads to this value being 3.
        Default [min,max]:  [0.1,4]
    
    See also
    --------
    crease_ga.Model.load_shape
    '''
    
    def __init__(self,
                 shape_params = [0.5,1],
                minvalu = (20, 0.1, 0.01, 0.01, -1, -1, -1, 0.1),
                maxvalu = (3000, 0.95, 0.99, 0.99, 1, 1, 1, 4)):
        scat_density = shape_params[0]
        sigmabead = shape_params[1]
        self._numvars = 2
        self.minvalu = minvalu
        self.maxvalu = maxvalu
        self.scat_density = scat_density    ## scatterer density
        self.sigmabead = sigmabead    ## scatterer density
    
    @property
    def numvars(self):
        return self._numvars
    

    def converttoIQ(self, qrange, param):
        '''
        Calculate computed scattering intensity profile.

        Parameters
        ----------
        qrange: numpy.array
            q values.
        param: numpy.array
            Decoded input parameters. See *Notes* section of the class
            documentation.

        Returns
        -------
        IQid: A numpy array holding I(q).
        '''
        # q values, decoded parameters, 
        # number of repeat units per chain, fraction of B beads per chain, core density, 
        # scatterer diameter, molar mass of B chemistry, 
        # length of A chemistry bond, length of B chemistry bond, 
        # number of scatterers per chain, # of replicates, stdev in Rcore size
        sigmabead = self.sigmabead
        scat_density = self.scat_density
        
        IQid=np.zeros((len(qrange)))      #initiates array for output IQ

        ### Parameters used to generate scatterer placements ###
        R_total=param[0]
        f_core=param[1]
        f_Ain=param[2]
        f_Aout=param[3]
        sld_Ain=param[4]     # split of type A scatterer 
        sld_B=param[5]     # split of type A scatterer 
        sld_Aout=param[6]   # variation in Rcore, dispersity
        #print(Rcore, dR_Ain, dR_B, dR_Aout, sAin)
        Background=10**(-param[6])

        R_core = R_total*f_core
        dR_Ain = R_total*(1-f_core)*f_Ain
        dR_Aout = f_Aout*(R_total-R_core-dR_Ain)
        dR_B = R_total-R_core-dR_Ain-dR_Aout
        
        r = genLP(R_core, dR_Ain, dR_B, dR_Aout, scat_density) 
        IQid = LPOmega(qrange, r, sigmabead, sld_Ain, sld_B, sld_Aout)
        maxIQ=np.max(IQid)                                  
        IQid=np.true_divide(IQid,maxIQ)                    # normalizes the I(q) to have its maximum = 1
        IQid+=Background                                   # add background
        return IQid
