import numpy as np
import numexpr as ne
def gen_scat_coords_cyl(nscat,length,diam):
    added = 0
    ret = []
    while added < nscat:
        x = np.random.random()*2-1
        y = np.random.random()*2-1
        if x**2+y**2 < 1:
            added += 1
            ret.append([x*diam/2,y*diam/2,np.random.random()*length])
    return np.array(ret)

def spherical_pq(r,q):
    return 36*np.pi*(4/3*np.pi*r**3*(np.sin(q*r)-q*r*np.cos(q*r))/(q*r)**3)**2

def write_xyz(fn,coords):
    f = open(fn,'w')
    f.write("{:d}\n".format(len(coords)))
    f.write("\n")
    for i in coords:
        f.write("0\t{:.5f}\t{:.5f}\t{:.5f}\n".format(i[0],i[1],i[2]))
    f.close()

def sq(pos,qs,pbc = False, box_length = 0):
    nbins = len(qs)
    sq = np.zeros((nbins))
    natoms = len(pos)
    for i in range(natoms):
        all_disp = pos[i,:]-pos[(i+1):,:]
        if pbc:
            all_disp = all_disp-box_length*np.round(all_disp/box_length)

        rij = np.sqrt(np.sum(np.square(all_disp),axis=1))

        #numexpr
        qi = qs[np.newaxis,:]
        R = rij[:,np.newaxis]
        sq = sq+np.sum(ne.evaluate("2*sin(R*qi)/(R*qi)"),axis=0)

    return sq

def rot_matrix(u,theta):
    ux = u[0]
    uy = u[1]
    uz = u[2]
    R = np.zeros((3,3))
    R[0,0] = np.cos(theta)+ux**2*(1-np.cos(theta))
    R[0,1] = ux*uy*(1-np.cos(theta))-uz*np.sin(theta)
    R[0,2] = ux*uz*(1-np.cos(theta))+uy*np.sin(theta)
    R[1,0] = uy*ux*(1-np.cos(theta))+uz*np.sin(theta)
    R[1,1] = np.cos(theta)+uy**2*(1-np.cos(theta))
    R[1,2] = uy*uz*(1-np.cos(theta))-ux*np.sin(theta)
    R[2,0] = uz*ux*(1-np.cos(theta))-uy*np.sin(theta)
    R[2,1] = uz*uy*(1-np.cos(theta))+ux*np.sin(theta)
    R[2,2] = np.cos(theta)+uz**2*(1-np.cos(theta))
    
    return R

def gen_snake(sigma, L, unitL):
    #unit vectors defining cross section:
    cs_unit = np.array([[1,0,0],[0,1,0]])
    #unit vector defining primary axis:
    ax_unit = np.array([0,0,1])
    coords = []
    css = []
    axs = []
    curr = np.array([0,0,0])
    for i in range(int(L/unitL)):
        coords.append(curr)
        css.append(cs_unit)
        axs.append(ax_unit)
        #generate random angles for vector
        v_theta = np.arccos(2*np.random.random()-1)
        v_phi = 2*np.pi*np.random.random()
        #convert from spherical to cartesian
        rot_u = np.array([np.cos(v_phi)*np.sin(v_theta),
                          np.sin(v_phi)*np.sin(v_theta),
                          np.cos(v_theta)])
        rot_theta = np.absolute(np.random.normal(0,sigma))
        if rot_theta > np.pi:
            rot_theta = np.pi
        rot_mat = rot_matrix(rot_u,rot_theta)
        ax_unit = np.matmul(rot_mat,ax_unit)
        curr = curr+ax_unit*unitL
        cs_unit = np.array([np.matmul(rot_mat,c) for c in cs_unit])
    return coords,axs,css

def gen_snake_scat(nscat,coords,css,diam):
    N = len(coords)
    scat_coords = []
    for i in range(nscat):
        interval = int(np.floor(np.random.random()*(N-1)))
        t = np.random.random()
        axis_pos = coords[interval]*t+coords[interval+1]*(1-t)
        dx = np.random.random()*2-1
        dy = np.random.random()*2-1
        while dx**2+dy**2 >= 1:
            dx = np.random.random()*2-1
            dy = np.random.random()*2-1
        scat_coords.append(axis_pos+css[interval][0]*dx*diam/2+css[interval][1]*dy*diam/2)
    return np.array(scat_coords)    
def gen_scat_coords_flexcyl(nscat,length,diam,sigma):
    coords,axs,css = gen_snake(sigma,length,1)
    return gen_snake_scat(nscat,coords,css,diam)

def compute_error(wil,predict,observe,metric='squared_log_error',observe_uncertainty=None):
    if observe_uncertainty == None:
        observe_uncertainty = np.sqrt(observe)
    metrics = ['squared_log_error','abs_chi2']
    if not metric in metrics:
        raise Exception('Invalid metric')
    if metric == 'squared_log_error':
        return wil*(np.log(np.true_divide(observe,predict)))**2  #squared log error 
    if metric == 'abs_chi2':
        return wil*((observe-predict)/observe_uncertainty)**2

