import numpy as np
import random
import numexpr as ne
import sys

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

def LPOmega(qrange, nAin, nAout, nB, r):                # qvalues number_of_B number_of_A scatterer_coordinates
    Ntot=nAin+nB+nAout                                  # Total number of scatterers to loop through
    omegaarrt=np.zeros((1,len(qrange)))                 # initiating array
    
    omegaarr=np.zeros((1,len(qrange)))              # initiating array
    rur=r[0,:,:]# selects
    rur=rur.transpose()
    for i in range(Ntot-1):                         # loops through index and all further indexes to prevent double counting 
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
        inds_double_check=np.argwhere(np.isnan(vals))
        if len(inds_double_check)>0:
            print('nan error!')
        vals = ne.evaluate("sum((vals), axis=0)")   # adds together scatterer contributions for each q value
        omegaarr+=vals
    omegaarr=np.true_divide(2*omegaarr,Ntot)+1      # 1 accounts for the guarenteed overlap of same bead  # 2* accounts for double counting avoided to reduce computational expense by looping for all other pairs
    omegaarrt+=omegaarr                             # stores values between loops

    return omegaarrt

def visualize(r, Rcore, dR_Ain, dR_B, dR_Aout, sigmabead):
    import py3Dmol
    view = py3Dmol.view()
    
    for ri in r[0,:,:].transpose():
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
    
           

def genLP(Rcore, dR_Ain, dR_B, dR_Aout, sigmabead, nAin, nAout, nB):  
        # core radius, inner A layer thickness, B layer thickness, outer A layer thickness, 
        # bead diameter, # of inner A beads, # of outer A beads, # of B beads

        ntot = nAin+nB+nAout
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

    
class scatterer_generator:
    '''
    The wrapper class for vesicle shape. Default length unit: Angstrom.
    

    Notes
    -----
    **The following 7 shape-specific descriptors are to be specified by user (see
    *Attributes*) as 
    a list, in the precise order as listed, while calling `Model.load_shape`
    to load this shape:**


    num_scatterers: 
        Number of scatterers used to represent a chain. Default: 24
    N: 
        Number of monomers in a chain. Default: 54
    eta_B:
        Packing fraction of scatterers in B layer. Default: 0.5
    lmono_b:
        Diameter of a monomer of chemistry B. Default: 50.4 A
    lmono_a:
        Diameter of a monomer of chemistry A. Default: 50.4 A
    fb: 
        Fraction of monomers in chain that are of B type. fa = 1-fb. Default: 0.55
    nLP:
        Number of replicates for each individual. Default: 7



    **The following 7 parameters are to be predicted, in the precise order
    as listed, by GA:**
    
    R_core:
        Core radius. Default [min,max]: [50 A, 400 A]
    t_Ain:
        Thickness of inner A layer. Default [min,max]: [30 A, 200 A]
    t_B:
        Thickness of B layer. Default [min,max]: [30 A, 200 A]
    t_Aout:
        Thickness of outer A layer. Default [min,max]: [30 A, 200 A]
    sigma_Ain:
        Split of solvophilic scatterers between inner and outer layers. 
        Default [min,max]: [0.1, 0.45]
    sigma_R:
        Dispersity in vesicle size as implemented in the core radius.
        Default [min,max]: [0.0, 0.45]
    log10(bg):
        Negative log10 of background intensity. 
        E.g. an background intensity of 0.001 leads to this value being 3.
        Default [min,max]:  [0.1,4]
    
    See also
    --------
    crease_ga.Model.load_shape
    '''
    
    def __init__(self,
                 shape_params = [24,54,0.5,50.4,50.4,0.55,7],
                minvalu = (50, 30, 30, 30, 0.1, 0.0, 0.1),
                maxvalu = (400, 200, 200, 200, 0.45, 0.45, 4)):
        num_scatterers = shape_params[0]
        N = shape_params[1]
        rho_B = shape_params[2]
        lmono_a = shape_params[3]
        lmono_b= shape_params[4]
        fb = shape_params[5]
        nLP = shape_params[6]
        self._numvars = 7
        self.minvalu = minvalu
        self.maxvalu = maxvalu
        self.num_scatterers=num_scatterers    ## number of scatterers per chain
        self.N=N                 ## Number of beads on chain
        self.rho_B=rho_B            ## density/volume fraction of beads in B layer
        
        self.lmono_a=lmono_a         ## Angstrom 'monomer contour length'
        self.lmono_b=lmono_b         ## Angstrom 'monomer contour length'
        self.MB=np.pi/6*(self.lmono_b)**3 ## volume of B monomer
        self.sigmabead=np.true_divide(self.N*self.lmono_b,self.num_scatterers) ## scatterer bead diameter
        self.fb=fb              ## fraction of B type monomers in chain

        self.nLP=nLP                ## number of replicates
    
    @property
    def numvars(self):
        return self._numvars
    
    def calculateScattering(self,qrange,params,output_dir,n_cores):
        '''
        Determine each individual's computed scattering intensity profile.

        Parameters
        ----------
        qrange: numpy.array
            q values.
        params: numpy.array
            Decoded input parameters. See *Notes* section of the class
            documentation.
        output_dir: string
            folder to output any generated files
        n_cores: int
            number of cores to use for parallelization (if implemented)

        Returns
        -------
        IQids: A numpy array holding each individual's I(q).
        '''
        IQids = []
        for val in range(len(params)):
            sys.stdout.write("\rindividual {:d}/{:d}".format(val+1,len(params)))
            sys.stdout.flush()
            IQid=self.converttoIQ(qrange, params[val])
            IQids.append(IQid)
        IQids = np.array(IQids)
        return IQids

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
        N = self.N
        fb = self.fb
        rho_B = self.rho_B
        MB = self.MB
        lmono_a = self.lmono_a
        lmono_b = self.lmono_b
        num_scatterers = self.num_scatterers
        nLP = self.nLP
        
        IQid=np.zeros((len(qrange)))      #initiates array for output IQ

        ### Parameters used to generate scatterer placements ###
        Rcore=param[0]
        dR_Ain=param[1]
        dR_B=param[2]
        dR_Aout=param[3]
        sAin=param[4]     # split of type A scatterer 
        sigmaR=param[5]   # variation in Rcore, dispersity
        #print(Rcore, dR_Ain, dR_B, dR_Aout, sAin)
        Background=10**(-param[6]) 

        varR = Rcore*sigmaR # variation in Rcore
        disper = np.array([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]) # fixed intervals of sigma
        sum_omegaarr=np.zeros((1,len(qrange)))

        for step in range(0, nLP):
            Rcore = param[0] + varR*disper[step + int((9-nLP)/2.)]          ## add displacement to Rcore
    #         print("disper = ", disper[step + int((9-nLP)/2.)])
    #         print("Rcore = ", Rcore)
            vol_B = (4/3.0)*np.pi*(np.power(Rcore + dR_Ain + dR_B, 3) 
                                   - np.power(Rcore + dR_Ain, 3)) ## volume of solvophobic layer B
            nagg = int(np.true_divide( rho_B*vol_B, N*fb*MB ))    ## number of chains in vesicle
            ntot = nagg*num_scatterers                            ## total number of scatterers
            nB = int(ntot*fb)                                     ## number of scatterers in B
            nAin = int(ntot*(1-fb)*sAin)                          ## number of scatterers in A_in
            nAout = int(ntot*(1-fb)*(1-sAin))                     ## number of scatterers in A_out

            for reps in range(0, 3):
                ### Generates scatterer positions in structure ###
                r = genLP(Rcore, dR_Ain, dR_B, dR_Aout, sigmabead, nAin, nAout, nB)

                ### Calculates omega from scatterers in shape ###
                sum_omegaarr += LPOmega(qrange, nAin, nAout, nB, r)

        omegaarr=np.true_divide(sum_omegaarr,nLP*3)        # average omega
        omegaarr=omegaarr.reshape(len(qrange),)
        Fqb=LPFbead(qrange,sigmabead)                      # calcualtes sphere shape factor
        F2qb=np.multiply(Fqb,Fqb)                          # Sphere shape factor square
        sqmm=np.ones((np.shape(Fqb)))                      # assuming dilute mixture the micelle-micelle structure factor = 1
        F2qb_sqmm=np.multiply(F2qb,sqmm)                   # determines the micelle form factor
        IQid=np.multiply(omegaarr,F2qb_sqmm)               # calculates Icomp
        maxIQ=np.max(IQid)                                  
        IQid=np.true_divide(IQid,maxIQ)                    # normalizes the I(q) to have its maximum = 1
        IQid+=Background                                   # add background
        return IQid
