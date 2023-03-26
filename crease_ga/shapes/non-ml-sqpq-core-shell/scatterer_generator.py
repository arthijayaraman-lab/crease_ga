import numpy as np
import numexpr as ne
from scipy import spatial, stats
#from crease_ga.exceptions import CgaError
from crease_ga import utils
from os import path
import multiprocessing as mp
from tensorflow import keras
import sys
from subprocess import check_output, run
from itertools import repeat

def generateScatterer(D, cc, rho):
    # ensure some randomness
    seed = int(np.random.random()*10000/(D*cc))
    np.random.seed(seed)
    nscatter = int(np.rint(rho * 4/3.*np.pi*(D/2)**3))
    nA = int(np.rint(rho * 4/3.*np.pi*(D/2*cc)**3))
    nB = nscatter -nA
    atype = np.ones((nscatter),dtype='int')
    atype[nA:]+= 1
    # random locations
    pos = np.random.random((nscatter,3))* 2 - 1
    norm = np.linalg.norm(pos,axis=1)
    pos = pos/norm[:,None]
    # give A and B scatterers different magnitude of displacement
    maga = np.random.random(nA) * D/2*cc
    magb = np.random.random(nB) * D/2*(1-cc) + D/2*cc
    pos[:nA,:] = pos[:nA,:] * maga[:,None]
    pos[nA:,:] = pos[nA:,:] * magb[:,None]
    atype = np.ones((nscatter),dtype='int')
    atype[nA:]+= 1
    return atype, pos

def micelleIqCalc(atype, pos, qrange, ffactor1, ffactor2):
    natoms = len(pos)
    nbins = len(qrange)
    I = np.zeros((nbins))
    Icore = np.zeros((nbins))
    Ishell = np.zeros((nbins))
    atype = atype.astype('int')
    # for atoms find the fi*fj*sin(qr)/(qr) for each pair
    ff = np.where(atype<2,ffactor1,ffactor2)
    for i in range(natoms):
        # pair distances
        rij = np.sqrt(np.sum(np.square(pos[i,:]-pos[(i+1):,:]),axis=1))
        # numexpr
        qi = qrange[np.newaxis, :]
        R = rij[:, np.newaxis]
        # factor of 2 necessary since only doing each pair once
        sq = ne.evaluate("2*sin(qi*R)/(qi*R)")
        a1 = int(atype[i]-1)
        # do the f^2 term in front
        if a1 < 1:
            fi = ffactor1
            Icore += fi*fi
            mask = (atype[(i+1):]-1 < 1)
            Icore += np.sum(fi*fi*sq[mask], axis=0)
            I += fi*fi
            I += np.sum(fi*ff[(i+1):][:,np.newaxis]*sq, axis=0)
        else:
            fi = ffactor2
            Ishell += fi*fi
            mask = (atype[(i+1):]-1 > 0)
            Ishell += np.sum(fi*fi*sq[mask], axis=0)
            I += fi*fi
            I += np.sum(fi*ff[(i+1):][:,np.newaxis]*sq, axis=0)
    # normalize I(q) to 1.0 at lowest q value
    I /= I[0]
    Icore /= Icore[0]
    Ishell /= Ishell[0]
    return I, Icore, Ishell

def adjustDomains(types, tree, pos, domain, cutoff, neighbors, id, dict):
    # method that deals with returning neighbors and adjusting domains
    # get all neighbors of all particles in domain
    neigh = np.array(tree.query_ball_point(pos[id], cutoff), dtype=int)
    # ignore any that have same domain number
    neigh = neigh[types[neigh] != domain]
    # need to check & convert any with other domain (and rerun query)
    change = np.unique(types[neigh])
    for i in range(len(change)-1):
        # get all of other domain's type1 neighbors and keep unique ones
        egg = np.array(dict[change[i+1]], dtype=int)
        neigh = np.append(neigh, egg[types[egg] == 0])
        # update types and dictionary
        del dict[change[i+1]]
        types[types == change[i+1]] = domain
    # add new neigh to neighbor list
    neighbors = np.append(neighbors, neigh.astype(int))
    # get unique ones
    neighbors = np.unique(neighbors).astype(int)
    neighbors = neighbors[types[neighbors] == 0]
    return types, neighbors

# generate particle distribution to match input parameters
def generateDistribution(phi1, D1, S1, D2, S2, natoms):
    nomN1 = int(np.rint(phi1*natoms*(D2/D1)**3/(1-phi1+phi1*(D2/D1)**3)))
    desVol1 = nomN1*4/3.*np.pi*(D1/2)**3
    binmin1 = stats.lognorm.ppf(0.01, S1, scale=D1)
    binmax1 = stats.lognorm.ppf(0.99, S1, scale=D1)
    sizes1 = np.linspace(binmin1, binmax1, 11)
    bins1 = np.append(sizes1-(sizes1[1]-sizes1[0])/2,
                      sizes1[-1]+(sizes1[1]-sizes1[0])/2)
    r = stats.lognorm.rvs(s=S1, loc=0, scale=D1, size=nomN1)
    it = 1
    while np.sum(((r < binmin1) | (r > binmax1))) > 0:
        it += 1
        r = np.where(((r < binmin1) | (r > binmax1)),
                     stats.lognorm.rvs(S1, size=1, scale=D1), r)
    hist1, edges = np.histogram(r, bins=bins1)
    tempVol1 = np.sum(hist1*(4.0/3.0)*np.pi*(sizes1/2.0)**3.0)
    err = abs(tempVol1-desVol1)
    errold = 0
    while err > 0.0001*desVol1:
        hist1 = hist1*desVol1/tempVol1
        hist1 = hist1.astype(int)
        tempVol1 = np.sum(hist1*(4.0/3.0)*np.pi*(sizes1/2.0)**3.0)
        errold = err
        err = abs(tempVol1-desVol1)
        if err == errold:
            break
    nomN2 = natoms - nomN1
    desVol2 = nomN2*4/3.*np.pi*(D2/2)**3
    binmin2 = stats.lognorm.ppf(0.01, S2, scale=D2)
    binmax2 = stats.lognorm.ppf(0.99, S2, scale=D2)
    sizes2 = np.linspace(binmin2, binmax2, 11)
    bins2 = np.append(sizes2-(sizes2[1]-sizes2[0])/2,
                      sizes2[-1]+(sizes2[1]-sizes2[0])/2)
    r = stats.lognorm.rvs(s=S2, loc=0, scale=D2, size=nomN2)
    it = 1
    while np.sum(((r < binmin2) | (r > binmax2))) > 0:
        it += 1
        r = np.where(((r < binmin2) | (r > binmax2)),
                     stats.lognorm.rvs(S2, size=1, scale=D2), r)
    hist2, edges = np.histogram(r, bins=bins2)
    tempVol2 = np.sum(hist2*(4.0/3.0)*np.pi*(sizes2/2.0)**3.0)
    err = abs(tempVol2-desVol2)
    errold = 0
    while err > 0.0001*desVol2:
        hist2 = hist2*desVol2/tempVol2
        hist2 = hist2.astype(int)
        tempVol2 = np.sum(hist2*(4.0/3.0)*np.pi*(sizes2/2.0)**3.0)
        errold = err
        err = abs(tempVol2-desVol2)
        if err == errold:
            break
    Nps = np.append(hist1, hist2)
    return Nps

def ps(self, param, individual,output_dir):
    D1 = param[0]
    S1 = param[1]
    D2 = param[2]
    S2 = param[3]
    phi1 = param[5]
    sq1 = param[6]
    sq2 = param[7]
    sq3 = param[8]
    # get diameter
    binmin1 = stats.lognorm.ppf(0.01, S1, scale=D1)
    binmax1 = stats.lognorm.ppf(0.99, S1, scale=D1)
    sizes1 = np.linspace(binmin1, binmax1, 11)
    binmin2 = stats.lognorm.ppf(0.01, S2, scale=D2)
    binmax2 = stats.lognorm.ppf(0.99, S2, scale=D2)
    sizes2 = np.linspace(binmin2, binmax2, 11)
    diameters = np.append(sizes1,sizes2)
    # set randomness
    seed = int(np.random.random()*10000/(1+individual))
    np.random.seed(seed)
    lens = generateDistribution(
        phi1, D1, S1, D2, S2, self.N)
    lens = np.array(lens)
    # generate grid pattern
    psize = np.max(diameters)
    pos = np.zeros((np.sum(lens), 3))
    atype = np.zeros((np.sum(lens)))
    side = int(np.sum(lens)**(1/3))
    k = 0
    flag = True
    while flag:
        # combined double loop
        for i, j in ((i1, j1) for i1 in range(side) for j1 in range(side)):
            spot1 = j + i*side + k*side*side
            pos[spot1, :] = [j*psize, i*psize, k*psize]
            # check if placed everything
            if spot1 + 1 == np.sum(lens):
                flag = False
                break
        k += 1
    # recenter positions
    center = np.average(pos, axis=0)
    pos -= center
    dict = {}
    tree = spatial.cKDTree(pos)
    # hold domain number, 0 is type A [unchanged]
    domain = 1
    types = np.zeros((len(atype)), dtype=int)
    # how many particles need to be type B during loop
    swap = np.sum(lens[:11])-1
    print('swaping', swap, 'out of', np.sum(lens))
    # randomly select a particle to be first domain
    old = np.random.choice(np.arange(len(atype)), size=1, replace=False)[0]
    types[old] = domain
    # initialize neighbor list for domains
    neighbors = []
    # produce first domain neighbors and update list
    cutoff = psize * 1.35
    new_p = tree.query_ball_point(pos[old], cutoff)
    # for each neighbor
    neighbors = []
    for j in range(len(new_p)):
        # ignore the particle previously removed
        if new_p[j] == old:
            continue
        # otherwise add to neigh list
        neighbors.append(new_p[j])
    neighbors = np.array(neighbors)
    dict[domain] = neighbors
    # now apply genes and go through rest of particles to swap
    for i in range(swap):
        # see if param[5] passes -- do we continue the current domain?
        if (np.random.random() <= sq1) and (len(neighbors) > 0):
            # pick a type1 neighbor of current domain to become type2
            while True:
                pick = np.random.choice(
                    np.arange(len(neighbors)), size=1, replace=False)[0]
                if types[neighbors[pick]] == 0:
                    # now check if want to weigh neighbor with more A/B contacts
                    neigh = tree.query_ball_point(pos[neighbors[pick]], cutoff)
                    bneigh = np.sum(types[neigh] > 0)
                    frac = bneigh/(len(neigh)-1)
                    if (np.random.random() <= 1-np.sqrt(np.square(sq2-frac))):
                        types[neighbors[pick]] = domain
                        break
                else:
                    # ignore any that have same domain number
                    print(types[neighbors[pick]], domain)
                    neighbors = neighbors[types[neighbors] == 0]
                    print('crap', len(neighbors), pick)
                    print(neighbors)
                    print(pos[types > 0], domain)

            # update neighbor list of domain and combine domains if needed
            types, neighbors = adjustDomains(
                types, tree, pos, domain, cutoff, neighbors, neighbors[pick], dict)
            neighbors = neighbors[types[neighbors] != domain]
            dict[domain] = neighbors
        else:
            # if param[5] fails, need to select spot to start next domain
            mask = types == 0
            locs = np.where(mask == True)[0]
            ranger = np.sum(mask)
            while True:
                pick = np.random.choice(
                    np.arange(ranger), size=1, replace=False)[0]
                # get neighbors to deal with param[7]
                neigh = tree.query_ball_point(pos[mask][pick], cutoff)
                bneigh = np.sum(types[neigh] > 0)
                frac = bneigh/(len(neigh)-1)
                # low gene2 means wants only SP, high means wants SMP
                if (np.random.random() <= 1-np.sqrt(np.square(sq3-frac))):
                    # update ntypes, domain
                    domain += 1
                    types[locs[pick]] = domain
                    # update neighbor list of domain and combine domains if needed
                    neighbors = np.array(tree.query_ball_point(
                        pos[locs[pick]], cutoff), dtype=int)

                    types, neighbors = adjustDomains(
                        types, tree, pos, domain, cutoff, neighbors, locs[pick], dict)
                    neighbors = neighbors[types[neighbors] != domain]
                    dict[domain] = neighbors
                    break
    uni, count = np.unique(types, return_counts=True)
    bb = np.ones((len(count)), dtype=bool)
    bb[0] = False
    count = count[bb]
    # update atype based on domains and randomly assigning diameters
    temp = np.sum(lens[:11])
    internal = np.arange(temp)
    # get random order of particle ids for type1
    c1 = np.random.choice(internal, size=temp, replace=False)
    # go through and assign first lens[i] to each type
    starter = 0
    ender = 0
    at1 = np.zeros((temp))
    for i in range(11):
        ender += lens[i]
        at1[c1[starter:ender]] = i+1
        starter = ender
    temp = np.sum(lens[11:])
    internal = np.arange(temp)
    # get random order of particle ids for type2
    c1 = np.random.choice(internal, size=temp, replace=False)
    starter = 0
    ender = 0
    # go through and assign first lens[i] to each type
    at2 = np.zeros((temp))
    for i in range(11):
        ender += lens[i+11]
        at2[c1[starter:ender]] = i+12
        starter = ender
    # finally update atype using at1, at2, and types
    atype[types == 0] = at2
    atype[types > 0] = at1
    # once all required particles are swapped to type2, need to generate lammps file
    v = 0
    # get particle volume
    for i in range(len(lens)):
        v += lens[i]*(diameters[i]/2.)**3
    eta = 0.25
    # move particles as much as possible
    length = np.max([np.max(pos), -np.min(pos)])
    # set length to get eta start of 0.58
    add = (4/3.*np.pi*v/eta)**(1/3.)
    l = add
    # here expand particle position to reduce overlap
    pos = pos + pos/(l/2.)*((l-length)/2.+(length/2.-np.max(pos)))
    f1 = open(output_dir+'datafile'+str(individual), 'w')
    header = [
        '#LAMMPS datafile\n{} atoms\n0 bonds\n22 atom types\n0 bond types\n0 angles\n'.format(len(atype))]
    header.append(
        '0 angle types\n0 dihedrals\n0 dihedral types\n0 impropers\n0 improper types\n\n')
    header.append('Masses\n\n')
    for i in range(22):
        header.append('{} 1.0\n'.format(i+1))
    header.append('\nAtoms\n\n')
    for i in range(len(header)):
        f1.write('{}'.format(header[i]))
    for i in range(len(atype)):
        f1.write('{} 0 {} {} {} {}\n'.format(
            i+1, int(atype[i]), pos[i, 0], pos[i, 1], pos[i, 2]))
    f1.close()
    le = np.abs(l/2.0)
    # need to write box size for lammps
    f1 = open(output_dir+'lammps.relax'+str(individual),'w')
    f1.write('units lj\natom_style bond\nboundary p p p\n')
    f1.write('region box block {} {} {} {} {} {} units box\n'.format(-le,le,-le,le,-le,le))
    f1.write('create_box 22 box\n')
    f1.write('print $(step) file check{}\npair_style lj/cut/opt {}\n'.format(individual,psize*2.5))
    for i in range(len(diameters)):
        for j in range(i,len(diameters)):
            temp = (diameters[i]+diameters[j])/2.0+20
            f1.write('pair_coeff {} {} 1.0 {} \n'.format(int(i+1),int(j+1),temp,temp*2**(1/6.)))
    f1.write('read_data datafile{} add append\nminimize 1.0e-8 1.0e-9 1000000 1000000\n'.format(individual))
    #f1.write('neighbor 66.0 bin\nrun_style verlet\n')
    f1.write('neighbor 66.0 bin\nrun_style verlet\ndump 1 all custom 5000 {}.txt id mol type x y z\ndump_modify 1 sort id append no flush yes\n'.format(individual))
    f1.write('pair_style colloid {}\n'.format(psize*2.5))
    for i in range(len(diameters)):
        for j in range(i,len(diameters)):
            f1.write('pair_coeff {} {} 0.1 1.0 {} {} {}\n'.format(int(i+1),int(j+1),diameters[i],diameters[j],(diameters[i]+diameters[j])/2*2.5))
    f1.write('variable xp atom -x/3000\nvariable yp atom -y/3000\nvariable zp atom -z/3000\n'.format(psize*2))
    f1.write('variable e atom x*x/6000+y*y/6000+z*z/6000\nfix move all addforce v_xp v_yp v_zp energy v_e\nfix_modify move energy yes\n')
    f1.write('minimize 1e-6 1.0e-7 1000000 1000000\nunfix move\nvariable xp atom -x/600\nvariable yp atom -y/600\nvariable zp atom -z/600\n')
    f1.write('fix move all addforce v_xp v_yp v_zp\nvelocity all create 1.0 {}\ntimestep 0.004\nfix tem all nvt temp 1.0 1.0 $(100*dt)\n'.format(np.random.randint(1,99999)))
    f1.write('print $(step) file check{}\ndump 2 all custom 5000 end{} id mol type x y z\ndump_modify 2 sort id append no flush yes\n run 30000\n'.format(individual,individual))
    f1.write('variable xp atom -x/300\nvariable yp atom -y/300\nvariable zp atom -z/300\n')
    f1.write('fix move all addforce v_xp v_yp v_zp\ntimestep 0.004\n'.format(np.random.randint(1,99999)))
    f1.write('print $(step) file check{}\ndump 2 all custom 3000 end{} id mol type x y z\ndump_modify 2 sort id append no flush yes\n run 42000\n'.format(individual,individual))
    f1.close()

def read(pair, d,output_dir):
    d = np.array(d)
    check = int(check_output(["bash",output_dir+"run3.sh",str(pair)]))
    if check == 0:
        return 0, 0, True
        
    with open(output_dir+'end'+str(pair), 'r') as f:
        text = f.read()
    lines = text.splitlines()
    # deal with failed individual
    if len(lines) < 2:
        return 0, 0, True
    if lines[0].split()[0] == "ERROR":
        return 0, 0, True
    n = int(lines[3])
    c = -1
    while True:
        t = lines[c].split()
        # check if atoms are same as start
        if t[-1] != 'ATOMS':
            c -= 1
            continue
        if int(lines[c+1]) != n:
            c -= 1
            continue
        # set starting line for below
        start = c + 7
        break
    pos = np.zeros((n, 3))
    atype = np.zeros((n), dtype=int)
    for i in range(n):
        line = lines[start+i].split()
        atype[i] = int(line[2])
        pos[i, 0] = float(line[3])
        pos[i, 1] = float(line[4])
        pos[i, 2] = float(line[5])
    # need to get outer particles -- voronoi
    r_outer = np.max(np.sum(np.square(pos), axis=1))
    vor = spatial.Voronoi(pos)
    vertices = vor.vertices
    regions = vor.regions
    region_id = vor.point_region
    # accept particles with vertices beyond furthest particle
    check = np.zeros(len(vertices))
    for i in range(len(vertices)):
        r = vertices[i, 0]**2+vertices[i, 1]**2+vertices[i, 2]**2
        if r > r_outer:
            check[i] = -1
    # now get ids of furthest particles
    id = []
    for i in range(len(region_id)):
        t = region_id[i]
        ee = regions[t]
        for j in range(len(ee)):
            if ee[j] == -1 or check[ee[j]] == -1:
                id.append(i)
                break
    poso = pos[id, :]
    ao = atype[id]
    # move to 1 radius * norm to get outer part of particle
    norm = np.linalg.norm(poso, axis=1)
    diam1 = np.array(d[ao-1], dtype=float)/2.
    temp = np.copy(poso)
    temp[:, 0] *= diam1/norm
    temp[:, 1] *= diam1/norm
    temp[:, 2] *= diam1/norm
    poso += temp
    vol = spatial.ConvexHull(poso).volume
    v = 0
    for i in range(len(d)):
        v += 4/3.*np.pi*len(atype[atype == (i+1)])*(d[i]/2)**3
    eta = v/vol
    # want eta > 0.5
    etacut = 0.5
    return pos, atype, eta < etacut

def setFormFactor(diameters, qrange, cs, sldc, slds,rho):
    fqtot = np.zeros((22,len(qrange)))
    fqcore = np.zeros((22,len(qrange)))
    fqshell = np.zeros((22,len(qrange)))
    for i in range(11):
        a,p = generateScatterer(diameters[i], cs, rho)
        fqtot[i], fqcore[i], fqshell[i] = micelleIqCalc(a, p, qrange, sldc,slds)
    return fqtot, fqcore, fqshell, np.zeros((11,len(qrange)))


def iqCalc(atype, pos, qrange, ffactor1, ffcore, ffshell, ffactor2):
    natoms = len(pos)
    nbins = len(qrange)
    I2 = np.zeros((nbins))
    I2core = np.zeros((nbins))
    I2shell = np.zeros((nbins))
    S2 = np.zeros((nbins))
    # for atoms find the fi*fj*sin(qr)/(qr) for each pair
    counter = 0
    for i in range(natoms):
        a1 = int(atype[i]-1)
        # do the f^2 term in front
        if a1 >10:
            counter += 1
            ### ignore type B particles as solvent
            continue
        else:
            fi = ffactor1[a1]
            I2 += fi*fi
            I2core += ffcore[a1]*ffcore[a1]
            I2shell += ffshell[a1]*ffshell[a1]
        # pair distances
        mask = (atype[(i+1):]-1 < 11)
        rij = np.sqrt(np.sum(np.square(pos[i,:]-pos[(i+1):,:][mask]),axis=1))
        #rij = np.sqrt(np.sum(np.square(pos[i,:]-pos[(i+1):,:]),axis=1))
        # numexpr
        qi = qrange[np.newaxis, :]
        R = rij[:, np.newaxis]
        # factor of 2 necessary since only doing each pair once
        sq = ne.evaluate("2*sin(qi*R)/(qi*R)")
        # apply form factor to get Iq
        # end with type B particles
        I2 += np.sum(ffactor1[a1]*ffactor1[atype[(i+1):][mask]-1]*sq, axis=0)
        I2core += np.sum(ffcore[a1]*ffcore[atype[(i+1):][mask]-1]*sq, axis=0)
        I2shell += np.sum(ffshell[a1]*ffshell[atype[(i+1):][mask]-1]*sq, axis=0)
        S2 += np.sum(sq,axis=0)
    S2 /= (natoms-counter)
    S2 += 1
    return I2, I2core, I2shell, S2

## this is for slit height = 0 & slit width = value
def weightMatrix(qext,q,h):
    # using SasView as template 
    qedge = np.hstack([qext[0]-0.5*(qext[1]-qext[0]),0.5*(qext[1:]+qext[:-1]),qext[-1]+0.5*(qext[-1]-qext[-2]),])
    weight = np.zeros((len(q),len(qext)))
    for i, qi in enumerate(q):
        # assumes w = 0
        inx = 1. * ((qext >= qi-h) & (qext <= qi+h))
        absx = 1. * (qext < abs(qi-h)) if qi<h else 0.
        weight[i,:] = (inx + absx) * np.diff(qedge) / (2.*h)
    return weight.transpose()

## this is for slit height = value & slit width = 0
def weightMatrix2(qext,q,h):
    qedge = np.hstack([qext[0]-0.5*(qext[1]-qext[0]),0.5*(qext[1:]+qext[:-1]),qext[-1]+0.5*(qext[-1]-qext[-2]),])
    weight = np.zeros((len(q),len(qext)))
    for i, qi in enumerate(q):
        u_limit = np.sqrt(qi**2 + h**2)
        u_edges = qedge**2 - qi**2
        u_edges[qedge < abs(qi)] = 0.
        u_edges[qedge > u_limit] = u_limit**2 - qi**2
        weight[i,:] = np.diff(np.sqrt(u_edges))/h
    return weight.transpose()

def smearSlit(q,qext,iq,h=0.23969):
    #weight = weightMatrix(qext,q,h)
    weight = weightMatrix2(qext,q,h)
    iqsmear = np.dot(iq[None,:],weight).flatten()
    return iqsmear

class scatterer_generator:
    def __init__(self,
                 shape_params=[1, 1,20000,False,0,'',''],
                 minvalu=(40,0.05,50,0.01,0.05,0.01,0,0,0,10),
                 maxvalu=(60,0.20,50,0.01,0.95,0.50,1,1,1,20)):
        
        # allow input of scattering differences between A and B type polymers
        sld1 = shape_params[0]
        sld2 = shape_params[1]
        self.sld1 = sld1
        self.sld2 = sld2
        # number of particles, sets box size
        N = shape_params[2]
        self.slit = shape_params[3]
        self.slit_h = shape_params[4]
        self.mlModelsq = shape_params[5]
        self.mlModelfq = shape_params[6]
        self.smear = self.slit
        self.N = N  # Number of particles to use
        
        # scatterer density for form factor calculation
        # larger improves resolution but increases calculation time    
        self.rho = 0.24
        
        # genes are:
            # 1) average micelle diameter
            # 2) micelle lognormal dispersity
            # 3) solvent bead diameter -- for spacing micelles out
            # 4) solvent lognormal dispersity
            # 5) A:B core:shell split -- fraction of total micelle diameter
            # 6) micelle concentration (vol)
            # 7) gene related to aggregation of micelles
            # 8) gene ralted to sphericity of aggregate
            # 9) gene related to aggreagte spacing
            # 10) background scattering
        self._numvars = 10
        self.minvalu = minvalu
        self.maxvalu = maxvalu
            
        if path.isfile('current_sse.txt'):
            self.best_fit = np.genfromtxt('current_sse.txt')
        else:
            self.best_fit = np.inf

    @property
    def numvars(self):
        return self._numvars
    
    def produceStructure(self, param, individual):
        ps(self, param, individual,self.output_dir)

    def calculateScattering(self,model,qrange,params,output_dir,n_cores=1):
        '''
        Determine each individual's computed scattering intensity profile.

        Parameters
        ----------
        qrange: numpy.array
            q values.
        params: numpy.array
            Decoded input parameters. See *Notes* section of the class
            documentation.
        n_cores: int
            number of CPU cores to use during multiprocessing

        Returns
        -------
        IQids: A numpy array holding each individual's I(q).
        '''
        if path.isfile(output_dir+'current_sse.txt'):
            self.best_fit = np.genfromtxt(output_dir+'current_fit.txt')
        else:
            self.best_fit = np.inf
        self.qrange = qrange
        self.output_dir = output_dir
        self.params = params
        ## produce structures for LAMMPS
        #"""
        pool = mp.Pool(n_cores)
        one = pool.starmap(self.produceStructure, 
                zip(params,range(len(params))))
        pool.close()
        pool.join()
        print('finished making datafiles')
        # now run LAMMPS to create close-packed strucutres
        f1 = open(output_dir+'run.sh','w')
        f1.write('d={}\ncd $d\nfor f in $( seq 0 $1 ); do\n  lmp -in ./lammps.relax$f > log.lammps &\n  echo "made" > "./done$f.txt"\ndone\nwait\ncd ..'.format(output_dir))
        f1.close()
        f1 = open(output_dir+'run2.sh','w')
        f1.write('d={}\necho $(ls $d/done* | wc -l)\n'.format(output_dir))
        f1.close()
        f1 = open(output_dir+'run3.sh','w')
        f1.write('d={}\necho $(wc -l "$d/$1.txt" | head -c 1)\n'.format(output_dir))
        f1.close()
        
        run(['bash',output_dir+'run.sh',str(len(params)-1)])
        while (int(check_output(['bash',output_dir+'run2.sh'])) < len(params)):
            time.sleep(30)
        
        print('finished LAMMPS')
        pool = mp.Pool(n_cores)
        one = pool.starmap(self.doScattering, zip(
                range(len(params)), params))
        pool.close()
        pool.join()
        print('finished scattering')
        one = np.array(one,dtype=object)
        self.atype = one[:, 0]
        self.pos = one[:, 1]
        self.params = params
        iq = one[:, 2]
        IQids = np.array(iq)
        print('I(q) not normalized')
        self.iq = IQids
        return IQids
    
    def doScattering(self, individual, param):
        D1 = param[0]
        S1 = param[1]
        D2 = param[2]
        S2 = param[3]
        binmin1 = stats.lognorm.ppf(0.01, S1, scale=D1)
        binmax1 = stats.lognorm.ppf(0.99, S1, scale=D1)
        sizes1 = np.linspace(binmin1, binmax1, 11)
        binmin2 = stats.lognorm.ppf(0.01, S2, scale=D2)
        binmax2 = stats.lognorm.ppf(0.99, S2, scale=D2)
        sizes2 = np.linspace(binmin2, binmax2, 11)
        diameters = np.append(sizes1,sizes2)
        p, a, flag = read(individual, diameters,self.output_dir)
        if flag:        # if close-packed failed, return large SSE
            print('error',individual)
            return [np.zeros((10000)), np.zeros((10000, 3)), -self.IQin]

        # need to set form factor stuff
        ff1, ff1core, ff1shell, ff2 = setFormFactor(diameters, self.qrange, param[4], self.sld1, self.sld2, self.rho)
        # save micelle form factor
        eggg, count = np.unique(a,return_counts=True)
        tempo = np.sum(ff1[:11]*count[:11,np.newaxis],axis=0)
        np.savetxt('fq-core-shell{}.txt'.format(individual),np.c_[tempo])

        # do I(q) scattering
        icomptot, icompcore, icompshell, temps = iqCalc(a, p, self.qrange, ff1, ff1core, ff1shell, ff2)
        np.savetxt('sq-core-shell{}.txt'.format(individual),np.c_[temps])
        
        Background = 10**(-param[-1])
        icomptot += Background
        icompcore += Background
        icompshell += Background
        icomp = np.append(icomptot,icompcore)
        icomp = np.append(icomp,icompshell)
        return [a, p, icomp]
        
    def postprocess(self,model):
        '''
        Save struture of best overall individual

        Parameters
        ----------
        model: object
            model.py self reference
        '''
        for i in range(len(self.pos)):
            pos = self.pos[i]
            atype = self.atype[i]
            iq = self.iq[i].flatten()
            best_param = self.params[i].flatten()
            # get diameter
            D1 = best_param[0]
            S1 = best_param[1]
            D2 = best_param[2]
            S2 = best_param[3]
            binmin1 = stats.lognorm.ppf(0.01, S1, scale=D1)
            binmax1 = stats.lognorm.ppf(0.99, S1, scale=D1)
            sizes1 = np.linspace(binmin1, binmax1, 11)
            binmin2 = stats.lognorm.ppf(0.01, S2, scale=D2)
            binmax2 = stats.lognorm.ppf(0.99, S2, scale=D2)
            sizes2 = np.linspace(binmin2, binmax2, 11)
            diameters = np.append(sizes1,sizes2)
            np.savetxt(self.output_dir+'current_diameters{}.txt'.format(i),np.c_[diameters])
            f1 = open(self.output_dir+'structure{}.txt'.format(i),'w')  # write as LAMMPS file
            f1.write('ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{}\nITEM: BOX BOUNDS pp pp pp\n'.format(len(pos)))
            f1.write('-3e5 3e5\n-3e5 3e5\n-3e5 3e5\n')
            f1.write('ITEM: ATOMS id mol type x y z\n')
            for j in range(len(pos)):
                f1.write('{} 0 {} {} {} {}\n'.format(
                    j+1, int(atype[j]), pos[j][0], pos[j][1], pos[j][2]))
            f1.close()
            f1 = open(self.output_dir+'icomp{}.txt'.format(i), 'w') 
            for j in range(len(iq)):
                f1.write('{} {} \n'.format(
                    model.qrange[j], iq[j]))
            f1.close()
        #continue
            
        # save necessary information for structure
        if np.min(model.fit) < self.best_fit:
            self.best_fit = np.min(model.fit)
            val = np.where(model.fit == self.best_fit)[0][0]
            pos = self.pos[val]
            atype = self.atype[val]
            iq = self.iq[val].flatten()
            best_param = self.params[val].flatten()
            del self.atype, self.pos, self.params, self.iq
            np.savetxt(self.output_dir+'current_fit.txt',np.c_[self.best_fit])
            np.savetxt(self.output_dir+'current_genes.txt',np.c_[best_param])
            # get diameter
            D1 = best_param[0]
            S1 = best_param[1]
            D2 = best_param[2]
            S2 = best_param[3]
            binmin1 = stats.lognorm.ppf(0.01, S1, scale=D1)
            binmax1 = stats.lognorm.ppf(0.99, S1, scale=D1)
            sizes1 = np.linspace(binmin1, binmax1, 11)
            binmin2 = stats.lognorm.ppf(0.01, S2, scale=D2)
            binmax2 = stats.lognorm.ppf(0.99, S2, scale=D2)
            sizes2 = np.linspace(binmin2, binmax2, 11)
            diameters = np.append(sizes1,sizes2)
            np.savetxt(self.output_dir+'current_diameters.txt',np.c_[diameters])
            f1 = open(self.output_dir+'best_structure.txt','w')  # write as LAMMPS file
            f1.write('ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{}\nITEM: BOX BOUNDS pp pp pp\n'.format(len(pos)))
            f1.write('-3e5 3e5\n-3e5 3e5\n-3e5 3e5\n')
            f1.write('ITEM: ATOMS id mol type x y z\n')
            for j in range(len(pos)):
                f1.write('{} 0 {} {} {} {}\n'.format(
                    j+1, int(atype[j]), pos[j][0], pos[j][1], pos[j][2]))
            f1.close()
            f1 = open(self.output_dir+'best_Icomp.txt', 'w') 
            for j in range(len(iq)):
                f1.write('{} {} \n'.format(
                    model.qrange[j], iq[j]))
            f1.close()
