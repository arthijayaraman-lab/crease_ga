#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 01:08:19 2021

@author: zijiewu
"""
from math import gamma
import numpy as np

def schulz_indpoint(x,mean,p,eval_method='log'):
    if eval_method == 'normal':
        z = (1-p**2)/p**2
        return (z+1)**(z+1)*(x/mean)**z*np.exp(-(z+1)*x/mean)/(mean*gamma(z+1))
    else:
        z = (1-p**2)/p**2
        return np.exp((z+1)*np.log(z+1)+z*np.log(x/mean)-(z+1)*x/mean-np.log(mean*gamma(z+1)))
def schulz_fxn(mean,p):
    z = (1-p**2)/p**2
    if mean-p*mean*4 < 0:
        xs = np.linspace(0,mean+p*mean*4,500)
    else:
        xs = np.linspace(mean-p*mean*4,mean+p*mean*4,500)
    fxs = np.zeros(len(xs))
    for xi,x in enumerate(xs):
        fxs[xi] = (z+1)**(z+1)*(x/mean)**z*np.exp(-(z+1)*x/mean)/(mean*gamma(z+1))
    fxs = fxs/np.sum(fxs)
    return xs,fxs

def schulz_rv(mean,p):
    xs,fxs = schulz_fxn(mean,p)
    pacc = np.cumsum(fxs)
    ri = np.random.uniform()
    return xs[np.where(pacc>=ri)[0][0]]

if __name__ == "main":
#if True:
    p=0.4
    xs,fxs = schulz_fxn(200,p)
    fig=plt.figure()
    plt.plot(xs,fxs)
    rvs = []
    for i in range(350):
        rvs.append(schulz_rv(200,p))
    bins,binedges = np.histogram(rvs,bins=20)
    plt.plot((binedges[:-1]+binedges[1:])/2,bins/np.sum(bins)/(binedges[1]-binedges[0]))
#    ax = plt.gca()
#    ax.set_xscale('log')
#    ax.set_yscale('log')
#    ax.set_ylim((1e-9,1e-1))
