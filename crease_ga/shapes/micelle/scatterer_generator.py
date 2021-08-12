import numpy as np
import random
import numexpr as ne

def LPFbead(qrange, sigmabead):
    '''
    Compute the spherical form factor given a range of q values.
    
    Parameters:
    ----------
    qrange: numpy.array
        array of values in q-space to compute form factor for.
    sigmabead: float
        diameter of the sphere.
    
    Return:
    ----------
    Fqb: numpy.array
        array of values of the spherical form factors (F(q)) computed at q-points listed in qrange.
    '''
    
    R=np.true_divide(sigmabead,2)
    QR=np.multiply(qrange,R)
    Fqb=np.multiply(np.true_divide(np.sin(QR)-np.multiply(QR,np.cos(QR)),np.power(QR,3)),3)  

    return Fqb

def LPOmega(qrange,nB,nA,nLP,r):
    Ntot=nB+nA
    omegaarrt=np.zeros((1,len(qrange)))
    for step in range(nLP):	
        omegaarr=np.zeros((1,len(qrange)))       
        rur=r[step,:,:]          
        for i in range(Ntot-1):
            x = np.square(rur[0,i]-rur[0,(i+1):])
            y = np.square(rur[1,i]-rur[1,(i+1):])
            z = np.square(rur[2,i]-rur[2,(i+1):])
            rij = np.sqrt(np.sum([x,y,z],axis=0))
            rs = rij[:,np.newaxis]
            Q = qrange[np.newaxis,:]
            vals = ne.evaluate("sin(Q*rs)/(Q*rs)")
            inds=np.argwhere(np.isnan(vals))
            if len(inds)>0:
                for val in inds:
                    vals[val[0],val[1]]=1
            inds_double_check=np.argwhere(np.isnan(vals))
            if len(inds_double_check)>0:
                print('nan error!')
            vals = ne.evaluate("sum((vals), axis=0)")
            omegatemp=vals#np.sum(vals,axis=1)
            omegaarr+=omegatemp
        omegaarr=np.true_divide(2*omegaarr,Ntot)+1 #1 accounts for the guarenteed overlap of same bead  #2* accounts for double counting avoided by looping for all other pairs
        omegaarrt+=omegaarr      
    omegaarrt=np.true_divide(omegaarrt,nLP)
    omegaarrt=omegaarrt.reshape(len(qrange),)
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

def genLP(Rcore,Rmicelle,sigmabead,nB,nA,nLP,Nagg):  #core radius, micelle radius, sigma to nm, #B beads, #A beads, #micelles
    Nfa=int(round(np.true_divide(nA,Nagg)))
    ntot=nA+nB
    power=2
    r=np.zeros((nLP,3,ntot))

    prb=np.zeros(100**3)
    rb=np.zeros((3,100**3))
    
    rsep=np.multiply(Rcore*0.01,range(100))
    v_theta=np.multiply(np.pi*0.01,range(100))
    v_phi=np.multiply(2*np.pi*0.01,range(100))
    
    vtheta, vphi = np.meshgrid(v_theta,v_phi)
    rax=np.sin(vtheta)*np.cos(vphi)
    ray=np.sin(vtheta)*np.sin(vphi)
    raz=np.cos(vtheta) 
            
    rax=rax.flatten()
    ray=ray.flatten()
    raz=raz.flatten() 
    
    Vol=np.power(rsep,power)
    
    trax, trsep = np.meshgrid(rax,rsep)
    tray, tVol = np.meshgrid(ray,Vol)
    traz, trsep = np.meshgrid(raz,rsep)
    
    rrax=np.multiply(trsep,trax)
    rray=np.multiply(trsep,tray)
    rraz=np.multiply(trsep,traz)
    
    rstorex=rrax.flatten()
    rstorey=rray.flatten()
    rstorez=rraz.flatten()

    rb=np.array((rstorex,rstorey,rstorez))       
    prb=tVol.flatten()
    #prb=np.tile(Vol, len(v_theta)*len(v_phi))          
    prb=np.cumsum(prb)
    prb=np.true_divide(prb,prb[-1])
    
    for step in range(nLP):  #creates configuration for each replicate nLP
        #Create initial random configuration
        for i in range(nB):
            pval=random.random()
            ind=np.where((prb>pval))[0][0]
            r[step,:,i]=rb[:,ind]


    pra=np.zeros(100**3)
    ra=np.zeros((3,100**3))
    
    rsep=np.multiply((Rmicelle-Rcore)*0.01,range(100))+Rcore
    v_theta=np.multiply(np.pi*0.01,range(100))
    v_phi=np.multiply(2*np.pi*0.01,range(100))
    
    vtheta, vphi = np.meshgrid(v_theta,v_phi)
    rax=np.sin(vtheta)*np.cos(vphi)
    ray=np.sin(vtheta)*np.sin(vphi)
    raz=np.cos(vtheta) 
            
    rax=rax.flatten()
    ray=ray.flatten()
    raz=raz.flatten() 
    
    Vol=np.power(rsep,power)
    
    trax, trsep = np.meshgrid(rax,rsep)
    tray, tVol = np.meshgrid(ray,Vol)
    traz, trsep = np.meshgrid(raz,rsep)
    
    rrax=np.multiply(trsep,trax)
    rray=np.multiply(trsep,tray)
    rraz=np.multiply(trsep,traz)
    
    rstorex=rrax.flatten()
    rstorey=rray.flatten()
    rstorez=rraz.flatten()

    ra=np.array((rstorex,rstorey,rstorez))       
    pra=tVol.flatten()
    #prb=np.tile(Vol, len(v_theta)*len(v_phi))          
    pra=np.cumsum(pra)
    pra=np.true_divide(pra,pra[-1])
    
    for step in range(nLP):  #creates configuration for each replicate nLP
        #Create initial random configuration
        for i in range(nA):
            pval=random.random()
            ind=np.where((pra>pval))[0][0]
            r[step,:,i+nB]=ra[:,ind]
             
    return r

    
class scatterer_generator:
    def __init__(self,
                 shape_params = [8,24,0.5,0.5,50.4, 50.4],
                minvalu = (2, 0, 0),
                maxvalu = (60, 1, 5)):
        ## genes are aggregation_number, ecorona, background
        num_scatterers = shape_params[0]
        N = shape_params[1] #polymer length #A+#B
        fa = shape_params[2] #'input fraction A'
        rho_core = shape_params[3] #'Density of solvophobic block' 
        lmono_a = shape_params[4] # Angstrom 'monomer contour length'
        lmono_b= shape_params[5] # Angstrom 'monomer contour length' 
        self._numvars = 3
        self.minvalu = minvalu
        self.maxvalu = maxvalu
        self.num_scatterers=num_scatterers    ## number of scatterers per chain
        self.N=N                 ## Number of beads on chain
        self.fa=fa              ## input fraction of A type polymer
        self.rho_core=rho_core            ## density/volume fraction of beads in core layer
        self.lmono_a=lmono_a         ## Angstrom 'monomer contour length'
        self.lmono_b=lmono_b         ## Angstrom 'monomer contour length'
        self.MB=np.pi/6*(self.lmono_b)**3 ## volume of B monomer
        self.sigmabead=np.true_divide(self.N*self.lmono_b,self.num_scatterers) ## scatterer bead diameter
    
    @property
    def numvars(self):
        return self._numvars
    

    def converttoIQ(self, qrange, param): 
        N = self.N
        fa = self.fa
        sigmabead = self.sigmabead
        rho_core = self.rho_core
        MB = self.MB
        lmono_a = self.lmono_a
        lmono_b = self.lmono_b
        num_scatterers = self.num_scatterers
        nLP=3  #3 replicates for determining omega
        IQid=np.zeros((len(qrange))) #output IQ
        nagg=int(param[0])                                       #why +1?
        Ecorona=param[1] 
        Background=10**(-param[2]) 
        nB=int(np.round(num_scatterers*(1-fa)*nagg))#
        nA=int(np.round(num_scatterers*fa*nagg))#      #nA number of beads in corona (nAi is beads per chain)
                           #8 beads per chain
        Rcore=0.5*np.power(np.true_divide(6*nagg*N*(1-fa)*MB,np.pi*rho_core),np.true_divide(1,3))
	 
        Rmicelle=Ecorona*N*fa*lmono_a    #this is the solvophillic backbone + side chain length (maximum extent possible)

        r=genLP(Rcore,Rmicelle,sigmabead,nB,nA,nLP,nagg) #generates scatterer positions in micelle
        omegaarr=LPOmega(qrange,nB,nA,nLP,r)        #calculates points that scatter in shape omega
        Fqb=LPFbead(qrange,sigmabead) #Sphere shape factor
        F2qb=np.multiply(Fqb,Fqb)
        sqmm=np.ones((np.shape(Fqb))) #assuming dilute mixture  LPSMMQ(qrange,etamicelle,Rmicelleprime) #micelle micelle structure factor
        F2qb_sqmm=F2qb#np.multiply(F2qb,sqmm)

        IQid=np.multiply(omegaarr,F2qb_sqmm)+Background  #calculates I from F2 sqmm = Fbead Fbead omega sqmm
        maxIQ=np.max(IQid)
        IQid=np.true_divide(IQid,maxIQ)

        return IQid


