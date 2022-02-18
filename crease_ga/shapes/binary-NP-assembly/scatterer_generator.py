import numpy as np
import random
import numexpr as ne
from scipy import spatial, stats
from crease_ga_diameter.exceptions import CgaError
import sys
from crease_ga_diameter import utils
from subprocess import check_output
from os import path


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


def ps(self, param, individual):
    phi1 = param[3]
    D1 = param[5]
    S1 = param[6]
    D2 = param[7]
    S2 = param[8]
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
    swap = np.sum(lens[-11:])-1
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
        # see if param[0] passes -- do we continue the current domain?
        if (np.random.random() <= param[0]) and (len(neighbors) > 0):
            # pick a type1 neighbor of current domain to become type2
            while True:
                pick = np.random.choice(
                    np.arange(len(neighbors)), size=1, replace=False)[0]
                if types[neighbors[pick]] == 0:
                    # now check if want to weigh neighbor with more A/B contacts
                    neigh = tree.query_ball_point(pos[neighbors[pick]], cutoff)
                    bneigh = np.sum(types[neigh] > 0)
                    frac = bneigh/(len(neigh)-1)
                    if (np.random.random() <= 1-np.sqrt(np.square(param[2]-frac))):
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
            # if param[0] fails, need to select spot to start next domain
            mask = types == 0
            locs = np.where(mask == True)[0]
            ranger = np.sum(mask)
            while True:
                pick = np.random.choice(
                    np.arange(ranger), size=1, replace=False)[0]
                # get neighbors to deal with param[1]
                neigh = tree.query_ball_point(pos[mask][pick], cutoff)
                bneigh = np.sum(types[neigh] > 0)
                frac = bneigh/(len(neigh)-1)
                # low gene2 means wants only SP, high means wants SMP
                if (np.random.random() <= 1-np.sqrt(np.square(param[1]-frac))):
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
    internal = np.arange(swap+1)
    # get random order of particle ids for type2
    c1 = np.random.choice(internal, size=(swap+1), replace=False)
    starter = 0
    ender = 0
    # go through and assign first lens[i] to each type
    at2 = np.zeros((swap+1))
    for i in range(11):
        ender += lens[i+11]
        at2[c1[starter:ender]] = i+12
        starter = ender
    # finally update atype using at1, at2, and types
    atype[types == 0] = at1
    atype[types > 0] = at2
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
    f1 = open('datafile'+str(individual), 'w')
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
    f1 = open('lammps.relax'+str(individual),'w')
    f1.write('units lj\natom_style bond\nboundary p p p\n')
    f1.write('region box block {} {} {} {} {} {} units box\n'.format(-le,le,-le,le,-le,le))
    f1.write('create_box 22 box\n')
    f1.write('print $(step) file check{}\npair_style lj/cut/opt {}\n'.format(individual,psize*2.5))
    for i in range(len(diameters)):
        for j in range(i,len(diameters)):
            temp = (diameters[i]+diameters[j])/2.0+20
            f1.write('pair_coeff {} {} 1.0 {} \n'.format(int(i+1),int(j+1),temp,temp*2**(1/6.)))
    f1.write('read_data datafile{} add append\nminimize 1.0e-8 1.0e-9 1000000 1000000\n'.format(individual))
    f1.write('neighbor 66.0 bin\nrun_style verlet\ndump 1 all custom 5000 {}.txt id mol type x y z\ndump_modify 1 sort id append no flush yes\n'.format(individual))
    f1.write('pair_style colloid {}\n'.format(psize*2.5))
    for i in range(len(diameters)):
        for j in range(i,len(diameters)):
            f1.write('pair_coeff {} {} 0.1 1.0 {} {} {}\n'.format(int(i+1),int(j+1),diameters[i],diameters[j],(diameters[i]+diameters[j])/2*2.5))
    #f1.write('variable xp atom -x/3000\nvariable yp atom -y/3000\nvariable zp atom -z/3000\n'.format(psize*2))
    #f1.write('variable e atom x*x/6000+y*y/6000+z*z/6000\nfix move all addforce v_xp v_yp v_zp energy v_e\nfix_modify move energy yes\n')
    #f1.write('minimize 1e-6 1.0e-7 1000000 1000000\nunfix move\nvariable xp atom -x/600\nvariable yp atom -y/600\nvariable zp atom -z/600\n')
    #f1.write('fix move all addforce v_xp v_yp v_zp\nvelocity all create 1.0 {}\ntimestep 0.004\nfix tem all nvt temp 1.0 1.0 $(100*dt)\n'.format(np.random.randint(1,99999)))
    #f1.write('print $(step) file check{}\ndump 2 all custom 5000 end{} id mol type x y z\ndump_modify 2 sort id append no flush yes\n run 30000\n'.format(individual,individual))
    f1.write('variable xp atom -x/300\nvariable yp atom -y/300\nvariable zp atom -z/300\n'.format(psize*2))
    f1.write('variable e atom x*x/600+y*y/600+z*z/600\nfix move all addforce v_xp v_yp v_zp energy v_e\nfix_modify move energy yes\n')
    f1.write('minimize 1e-6 1.0e-7 1000000 1000000\nunfix move\nvariable xp atom -x/300\nvariable yp atom -y/300\nvariable zp atom -z/300\n')
    f1.write('fix move all addforce v_xp v_yp v_zp\nvelocity all create 1.0 {}\ntimestep 0.002\nfix tem all nvt temp 1.0 1.0 $(100*dt)\n'.format(np.random.randint(1,99999)))
    #f1.write('print $(step) file check{}\nrun 35000\ndump 2 all custom 1 best-pos2.txt id mol type x y z\nrun 0'.format(individual))
    f1.write('print $(step) file check{}\ndump 2 all custom 5000 end{} id mol type x y z\ndump_modify 2 sort id append no flush yes\n run 35000\n'.format(individual,individual))
#    f1.write('variable xp atom -x/300\nvariable yp atom -y/300\nvariable zp atom -z/300\n')
#    f1.write('fix move all addforce v_xp v_yp v_zp\ntimestep 0.0035\n'.format(np.random.randint(1,99999)))
#    f1.write('print $(step) file check{}\ndump 2 all custom 3000 end{} id mol type x y z\ndump_modify 2 sort id append no flush yes\n run 60000\n'.format(individual,individual))
    f1.close()


def read(gen, pair, d):
    d = np.array(d)
    check = int(check_output(["bash","run3.sh",str(pair)]))
    if check == 0:
        return 0, 0, True
        
    with open('end'+str(pair), 'r') as f:
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
    #print('Gen {} individual {}: eta = {}'.format(gen, pair, eta))
    # want eta > 0.5
    etacut = 0.54
    return pos, atype, eta < etacut


def setFormFactor(fFactor, diameters, qrange, cs_split, custom_form):
    builtin_forms = ["hard_sphere", "core_shell", "custom"]
    if not (fFactor in builtin_forms):
        print(fFactor, builtin_forms)
        print(fFactor in builtin_forms)
        raise CgaError(
            "Form factor setting is not defined. Choose another or 'custom' to directly set the form factor")
    if fFactor in builtin_forms[0]:
        rq = diameters[:, np.newaxis]/2. * qrange[np.newaxis, :]
        volumes = 4./3.*np.pi*(diameters/2.)**3
        f = 3.0*volumes[:, np.newaxis]*(np.sin(rq)-rq*np.cos(rq))/(rq*rq*rq)
        f1 = np.copy(f)
        f1[11:] = 0
        f2 = np.copy(f)
        f2[:11] = 0
        return f1, f2
    elif fFactor in builtin_forms[1]:
        # core-shell morphologies are A) type2-type1 and B) type1-type2
        rqOut = diameters[:, np.newaxis]/2. * qrange[np.newaxis, :]
        volOut = 4./3.*np.pi*(diameters/2.)**3
        fOut = 3.0*volOut[:, np.newaxis] * \
            (np.sin(rqOut)-rqOut*np.cos(rqOut))/(rqOut*rqOut*rqOut)
        rqIn = cs_split * diameters[:, np.newaxis]/2. * qrange[np.newaxis, :]
        volIn = 4./3.*np.pi*(cs_split * diameters/2.)**3
        fIn = 3.0*volIn[:, np.newaxis] * \
            (np.sin(rqIn)-rqIn*np.cos(rqIn))/(rqIn*rqIn*rqIn)
        f1 = np.copy(fOut)
        # f1 gets outer[:11] and inner [11:]
        f1[11:] = fIn[11:]
        f2 = np.copy(fOut)
        # f2 gets outer[11:] and inner [:11]
        f2[:11] = fIn[:11]
        return f1, f2
    elif fFactor in builtin_forms[2]:
        # custom assumes user is aware of proper format for form factor
        # to ensure correct scattering calculation
        loadvals = np.loadtxt(custom_form)
        qload = loadvals[:, 0]
        FQin1 = loadvals[:, 1]
        FQin2 = loadvals[:, 2]
        # interpolate to ensure qrange has form factor at each q[i]
        f1 = np.interp(qrange, qload, FQin1)
        f2 = np.interp(qrange, qload, FQin2)
        return f1, f2


def iqCalc(atype, pos, qrange, ffactor1, ffactor2):
    natoms = len(pos)
    nbins = len(qrange)
    I1 = np.zeros((nbins))
    I2 = np.zeros((nbins))
    # for atoms find the fi*fj*sin(qr)/(qr) for each pair
    for i in range(natoms):
        a1 = int(atype[i]-1)
        # do the f^2 term in front
        if a1 < 11:
            fi = ffactor1[a1]
            I1 += fi*fi
        else:
            fi = ffactor2[a1]
            I2 += fi*fi
        # pair distances
        rij = np.sqrt(np.sum(np.square(pos[i,:]-pos[(i+1):,:]),axis=1))
        # numexpr
        qi = qrange[np.newaxis, :]
        R = rij[:, np.newaxis]
        # factor of 2 necessary since only doing each pair once
        sq = ne.evaluate("2*sin(qi*R)/(qi*R)")
        # apply form factor to get Iq
        # start with type A particles
        mask = (atype[(i+1):]-1 < 11)
        I1 += np.sum(ffactor1[a1]*ffactor1[atype[(i+1):]-1][mask]*sq[mask], axis=0)
        # end with type B particles
        mask = (atype[(i+1):]-1 >= 11)
        I2 += np.sum(ffactor2[a1]*ffactor2[atype[(i+1):]-1][mask]*sq[mask], axis=0)
    # normalize I(q) to 1.0 at highest q value
    #I1 /= I1[-1]
    return I1, I2

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
                 chemistry_params=['hard_sphere', 20000,False,0],
                 minvalu=(0, 0, 0, 0, 0.001,200,0.,200,0.),
                 maxvalu=(1, 1, 1, 1, 4,240,0.15,240,0.15), core_shell_split=0.5, custom_form=None,pinholeSmear=False,slitSmear=False,slit_h=0):
        form_factor = chemistry_params[0]
        N = chemistry_params[1]
        self.slit = chemistry_params[2]
        self.slit_h = chemistry_params[3]
        self.smear = self.slit
        self._numvars = 9
        self.minvalu = minvalu
        self.maxvalu = maxvalu
        self.form_factor = form_factor  # form factor of particles
        self.N = N  # Number of particles to use
        self.cs_split = core_shell_split
        self.cust_form = custom_form
        if path.isfile('current_sse.txt'):
            self.best_sse = np.genfromtxt('current_sse.txt')

    @property
    def numvars(self):
        return self._numvars
    
    def setIQload(self,qrange,IQin,IQin2,IQerr,IQerr2,popnumber,generations,nloci):
        self.qrange = qrange
        self.IQin = IQin
        self.IQin2 = IQin2
        self.IQerr = IQerr
        self.IQerr2 = IQerr2
        self.popnumber = popnumber
        self.generations = generations
        self.nloci = nloci

    def produceStructure(self, param, individual):
        ps(self, param, individual)
        return 0

    def doScattering(self, pop, individual, param):
        Background = 10**(-param[4])
        D1 = param[5]
        S1 = param[6]
        D2 = param[7]
        S2 = param[8]
        binmin1 = stats.lognorm.ppf(0.01, S1, scale=D1)
        binmax1 = stats.lognorm.ppf(0.99, S1, scale=D1)
        sizes1 = np.linspace(binmin1, binmax1, 11)
        binmin2 = stats.lognorm.ppf(0.01, S2, scale=D2)
        binmax2 = stats.lognorm.ppf(0.99, S2, scale=D2)
        sizes2 = np.linspace(binmin2, binmax2, 11)
        diameters = np.append(sizes1,sizes2)
        p, a, flag = read(pop, individual, diameters)
        if flag:        # if close-packed failed, return large SSE
            print('error',individual)
            return [np.zeros((10000)), np.zeros((10000, 3)), -self.IQin,-self.IQin2]

        if self.smear == True:
            if self.slit == True:
                h = self.slit_h
                qmin, qmax = np.min(self.qrange-h), np.max(np.sqrt((self.qrange-h)**2))
                if qmin < 0:
                    qmin = self.qrange[0]*0.02
                log_delta_q = (len(self.qrange) - 1) / (np.log(self.qrange[-1]) - np.log(self.qrange[0]))
                nlow = int(np.ceil(log_delta_q * (np.log(self.qrange[0])-np.log(qmin))))
                qlow = np.logspace(np.log10(qmin), np.log10(self.qrange[0]), nlow+1)[:-1]
                nhigh = int(np.ceil(log_delta_q * (np.log(qmax)-np.log(self.qrange[-1]))))
                qhigh = np.logspace(np.log10(self.qrange[-1]), np.log10(qmax), nhigh+1)[1:]
                qext = np.concatenate([qlow,self.qrange,qhigh])
                ff1, ff2 = setFormFactor(
                    self.form_factor, diameters, qext, self.cs_split, self.cust_form)
                icomp = iqCalc(a, p, qext, ff1, ff2)
                icompSmear = smearSlit(self.qrange, qext, icomp, h)
                icompSmear = np.true_divide(icompSmear,icompSmear[0])
                icompSmear += Background
                return [a, p, icompSmear]
                

        else:
            # need to set form factor stuff
            ff1, ff2 = setFormFactor(
                self.form_factor, diameters, self.qrange, self.cs_split, self.cust_form)
            icomp,icomp2 = iqCalc(a, p, self.qrange, ff1, ff2)
            icomp = np.true_divide(icomp,icomp[0])
            icomp2 = np.true_divide(icomp2,icomp2[0])
            icomp += Background
            icomp2 += Background
            return [a, p, icomp,icomp2]

    def fitness(self, pop, generation, output_dir, atype, pos, iq,iq2, params, metric='chi2'):
        cs = 10
        F1 = open(output_dir+'z_temp_results_'+str(generation)+'.txt', 'w')
        F1.write('Individual gene1 gene2 gene3 composition background error\n')
        #np.savetxt(output_dir+'z_temp_population_' +
        #           str(generation)+'.txt', np.c_[pop])

        fitn = np.zeros(self.popnumber)
        fitnfr = np.zeros(self.popnumber)
        fit = np.zeros(self.popnumber)
        qfin = self.qrange[-1]
        IQid_str = []
        params = []
        for val in range(self.popnumber):
            # gets the current structure variables
            param = utils.decode(pop, val, self.nloci,
                                 self.minvalu, self.maxvalu)
            params.append(param)
            IQid = iq[val]
            IQid2 = iq2[val]

            err = 0
            for qi, qval in enumerate(self.qrange):
                if (IQid[qi] > 0) & (self.IQin[qi] > 0):
                    if metric == 'log_sse':
                        if (qi < qfin):
                            # weighting factor
                            wil = np.log(np.true_divide(
                                self.qrange[qi+1], self.qrange[qi]))
                        else:
                            # weighting factor
                            wil = np.log(np.true_divide(
                                self.qrange[qi], self.qrange[qi-1]))
                        # squared log error
                        err += wil * \
                            (np.log(np.true_divide(
                                self.IQin[qi], IQid[qi])))**2
                    elif metric == 'chi2':
                        if self.IQerr is None:
                            # chi^2 with weighting of sqrt(IQin)
                            err += np.true_divide(
                                np.square(self.IQin[qi]-IQid[qi]), np.square(self.IQin[qi]))
                            err += np.true_divide(
                                np.square(self.IQin2[qi]-IQid2[qi]), np.square(self.IQin2[qi]))
                        else:
                            # chi^2 with weighting of (IQerr)**2
                            err += np.true_divide(
                                np.square(self.IQin[qi]-IQid[qi]), np.square(self.IQerr[qi]))
                else:
                    if self.IQerr is None:
                        err += 2 #*self.IQin[qi]/np.sqrt(self.IQin[qi])
                    else:
                        err += 2*self.IQin[qi]/self.IQerr[qi]
            fit[val] = err
            IQid_str.append(IQid)
            F1.write(
                ('{:>2d} '.format(val)+'{:.4f} '.format(param[0])+'{:.4f} '.format(param[1])+'{:.4f} '.format(param[2])+'{:.4f} '.format(param[3])+'{:.4f} '.format(param[4])+'{:.4f} '.format(param[7])+'{:.4f} '.format(param[8])+'{:.5g} '.format(err)+'\n'))
            F1.flush()

            """
            print('outputing all information to check error calculation')
            f21 = open("xhelp-pos-"+str(val),'w')
            D1 = params[val][5]
            S1 = params[val][6]
            D2 = params[val][7]
            S2 = params[val][8]
            binmin1 = stats.lognorm.ppf(0.01, S1, scale=D1)
            binmax1 = stats.lognorm.ppf(0.99, S1, scale=D1)
            sizes1 = np.linspace(binmin1, binmax1, 11)
            binmin2 = stats.lognorm.ppf(0.01, S2, scale=D2)
            binmax2 = stats.lognorm.ppf(0.99, S2, scale=D2)
            sizes2 = np.linspace(binmin2, binmax2, 11)
            diameters = np.append(sizes1,sizes2)
            f2dia = open('xhelp-diameters-'+str(val),'w')
            for i in range(len(diameters)):
                f2dia.write('{}\n'.format(diameters[i]))
            f2dia.close()
            f21.write('ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{}\nITEM: BOX BOUNDS pp pp pp\n'.format(
                len(atype[val].flatten())))
            f21.write('-3e5 3e5\n-3e5 3e5\n-3e5 3e5\n')
            f21.write('ITEM: ATOMS id mol type x y z\n')
            for j in range(len(atype[val].flatten())):
                if len(pos[val].shape) == 3:
                    f21.write('{} 0 {} {} {} {}\n'.format(j+1, int(atype[val].flatten()[j]), pos[val][0][j][0], pos[val][0][j][1], pos[val][0][j][2]))
                else:
                    f21.write('{} 0 {} {} {} {}\n'.format(j+1, int(atype[val].flatten()[j]), pos[val][j][0], pos[val][j][1], pos[val][j][2]))
            f21.close()
            f21 = open('xhelp-icomp-'+str(val), 'w')  # I(q)
            for j in range(len(IQid)):
                f21.write('{} {} \n'.format(
                    self.qrange[j], iq[val].flatten()[j]))
            f21.close()
            """

    
        maxerr = np.max(fit)  # determines maximum SSerror for the population
        fitn = np.subtract(maxerr,fit) # determines error differences
        bestfit = np.max(fitn)
        sumup = np.sum(fitn)
        avgfit = np.true_divide(sumup, self.popnumber)
        dval = bestfit-avgfit
        # linear scaling with cs as a scaleFactor
        ascale = np.true_divide(avgfit, dval)*(cs-1.0)
        bscale = avgfit*(1.0-ascale)
        sumup = 0
        # get scaled fitness to enable selection of bad candidates
        for val in range(self.popnumber):
            if (fitn[val] > avgfit):
                fitnfr[val] = ascale*fitn[val]+bscale
            else:
                fitnfr[val] = fitn[val]
        sumup = np.sum(fitnfr)
        pacc = np.zeros(self.popnumber)
        prob = np.true_divide(fitnfr, sumup)
        pacc = np.cumsum(prob)
        ### returns cummulative relative error from which individuals can be selected ###
        maxfit = np.min(fit)
        elitei = np.where(fit == maxfit)[0]                  # Best candidate
        secondfit = sorted(fit)[1]
        # Second best candidate
        secondi = np.where(fit == secondfit)[0]
        avgfit = np.average(fit)
        avgi = np.array([(np.abs(fit-avgfit)).argmin()])   # Average candidate
        minfit = np.max(fit)
        mini = np.where(fit == minfit)[0]                    # Worst candidate
        if avgfit == 0:
            avgfit = 1
        gdm = np.true_divide(maxfit, avgfit)
        if len(elitei) > 1:
            elitei = elitei[0]
        if len(secondi) > 1:
            secondi = secondi[0]
        if len(avgi) > 1:
            avgi = avgi[0]
        if len(mini) > 1:
            mini = mini[0]
        # for structure, need to save atype, pos, & iq as the structure matters
        params = np.array(params)
        if generation == 0:
            self.best_pos = pos[elitei][0]
            self.best_atype = atype[elitei][0]
            self.best_iq = iq[elitei][0]
            self.best_iq2 = iq2[elitei][0]
            self.best_sse = fit[elitei]
            print('params',params.shape)
            print('elitei',elitei)
            print('both',params[elitei])
            best_param = params[elitei][0]
            np.savetxt('current_sse.txt',np.c_[fit[elitei]])
            np.savetxt('current_genes.txt',np.c_[best_param])
            # get diameter
            D1 = best_param[5]
            S1 = best_param[6]
            D2 = best_param[7]
            S2 = best_param[8]
            binmin1 = stats.lognorm.ppf(0.01, S1, scale=D1)
            binmax1 = stats.lognorm.ppf(0.99, S1, scale=D1)
            sizes1 = np.linspace(binmin1, binmax1, 11)
            binmin2 = stats.lognorm.ppf(0.01, S2, scale=D2)
            binmax2 = stats.lognorm.ppf(0.99, S2, scale=D2)
            sizes2 = np.linspace(binmin2, binmax2, 11)
            diameters = np.append(sizes1,sizes2)
            np.savetxt('current_diameters.txt',np.c_[diameters])
            f1 = open(output_dir+'best_structure.txt',
                      'w')  # write as LAMMPS file
            f1.write('ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{}\nITEM: BOX BOUNDS pp pp pp\n'.format(
                len(pos[elitei][0])))
            f1.write('-3e5 3e5\n-3e5 3e5\n-3e5 3e5\n')
            f1.write('ITEM: ATOMS id mol type x y z\n')
            print('elitei',elitei)
            print('pos',pos[elitei].shape)
            print('atype',atype[elitei].shape)
            print('iq',iq[elitei].shape)
            for j in range(len(pos[elitei][0])):
                f1.write('{} 0 {} {} {} {}\n'.format(
                    j+1, int(atype[elitei][0][j]), pos[elitei][0][j][0], pos[elitei][0][j][1], pos[elitei][0][j][2]))
            f1.close()
            f1 = open(output_dir+'best_Icomp.txt', 'w')  # I(q)
#            f1.write('q Icomp(q) \n')
            for j in range(len(IQid)):
                f1.write('{} {} {}\n'.format(
                    self.qrange[j], iq[elitei][0][j],iq2[elitei][0],[j]))
            f1.close()
        else:
            if fit[elitei] < self.best_sse:
                self.best_pos = pos[elitei]
                self.best_atype = atype[elitei]
                self.best_iq = iq[elitei]
                self.best_iq2 = iq2[elitei]
                self.best_sse = fit[elitei]
                best_param = params[elitei][0]
                np.savetxt('current_sse.txt',np.c_[fit[elitei]])
                np.savetxt('current_genes.txt',np.c_[best_param])
                # get diameter
                D1 = best_param[5]
                S1 = best_param[6]
                D2 = best_param[7]
                S2 = best_param[8]
                binmin1 = stats.lognorm.ppf(0.01, S1, scale=D1)
                binmax1 = stats.lognorm.ppf(0.99, S1, scale=D1)
                sizes1 = np.linspace(binmin1, binmax1, 11)
                binmin2 = stats.lognorm.ppf(0.01, S2, scale=D2)
                binmax2 = stats.lognorm.ppf(0.99, S2, scale=D2)
                sizes2 = np.linspace(binmin2, binmax2, 11)
                diameters = np.append(sizes1,sizes2)
                np.savetxt('current_diameters.txt',np.c_[diameters])
                f1 = open(output_dir+'best_structure.txt',
                          'w')  # write as LAMMPS file
                f1.write('ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{}\nITEM: BOX BOUNDS pp pp pp\n'.format(
                    len(pos[elitei][0])))
                f1.write('-3e5 3e5\n-3e5 3e5\n-3e5 3e5\n')
                f1.write('ITEM: ATOMS id mol type x y z\n')
                for j in range(len(pos[elitei][0])):
                    f1.write('{} 0 {} {} {} {}\n'.format(
                        j+1, int(atype[elitei][0][j]), pos[elitei][0][j][0], pos[elitei][0][j][1], pos[elitei][0][j][2]))
                f1.close()
                f1 = open(output_dir+'best_Icomp.txt', 'w')  # I(q)
#                f1.write('q Icomp(q)\n')
                for j in range(len(IQid)):
                    f1.write('{} {} {}\n'.format(
                        self.qrange[j], iq[elitei][0][j], iq2[elitei][0][j]))
                f1.close()

        f = open(output_dir+'fitness_vs_gen.txt', 'a')
        if generation == 0:
            f.write('gen mini min avgi avg secondi second besti best\n')
        f.write('%d ' % (generation))
        f.write('%d %.8lf ' % (mini, minfit))
        f.write('%d %.8lf ' % (avgi, avgfit))
        f.write('%d %.8lf ' % (secondi, secondfit))
        f.write('%d %.8lf ' % (elitei, maxfit))
        f.write('\n')
        f.close()
        print('Generation best fitness: {:.4f}'.format(maxfit))
        print('Generation best fitness: {:.3f}'.format(gdm))
        params = np.array(params)
        print('Generation best parameters '+str(params[elitei]))

        return pacc, gdm, elitei, IQid_str
