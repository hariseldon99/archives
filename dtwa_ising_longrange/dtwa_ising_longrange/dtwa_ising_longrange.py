#!/usr/bin/env python
#Author: Analabha Roy

from __future__ import division, print_function
from mpi4py import MPI
from reductions import Intracomm
from redirect_stdout import stdout_redirected

import sys
import copy
import random
import numpy as np
from scipy.signal import fftconvolve
from scipy.sparse import *

from scipy.integrate import odeint

from pprint import pprint
from tabulate import tabulate

threshold = 1e-4
root = 0
#This is the kronecker delta symbol for vector indices
deltaij = np.eye(3)

#This is the Levi-Civita symbol for vector indices
eijk = np.zeros((3, 3, 3))
eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] = 1
eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1

def t_deriv(quantities, times):
    dt = np.gradient(times)
    return np.gradient(quantities, dt)

def drive(t, params):
    return params.h0 * np.cos(params.omega * t)

def weyl_hamilt(s,times,param):
    """
    Evaluates the Weyl Symbols of the Hamiltonian, H_w
    Does this at all times
    If |s^a> = (s^a_0, s^a_1 ... s^a_N), and
    H_w = -(1/2) * \sum_{nm} J_{nm} (J_x s^n_x s^m_x + J_y s^n_y s^m_y
            + J_z s^n_z s^m_z) - h(t) * \sum_n (h_x s^n_x +h_y s^n_y
            + h_z s^n_z)
    """
    N = param.latsize
    #s[:, 0:N] = sx , s[:, N:2*N] = sy, s[:, 2*N:3*N] = sz
    drvs = drive(times, param)
    hw = param.jx * np.dot(s[:,0*N:1*N],param.jmat.dot(s[:,0*N:1*N].T))
    hw += param.jy * np.dot(s[:,1*N:2*N],param.jmat.dot(s[:,1*N:2*N].T))
    hw += param.jz * np.dot(s[:,2*N:3*N],param.jmat.dot(s[:,2*N:3*N].T))
    hw = hw /(2.0 * param.norm)
    hw += drvs * (param.hx * np.sum(s[:, 0:N]) +\
      param.hy * np.sum(s[:, N:2*N]) + param.hz * np.sum(s[:, 2*N:3*N]))
    return -hw

def func_1storder(s, t, param):
    """
    The RHS of general case, per Schachemmayer eq A2
    """
    N = param.latsize
    #s[0:N] = sx , s[N:2*N] = sy, s[2*N:3*N] = sz
    drv = drive(t, param)
    jsx = 2.0 * param.jx * param.jmat.dot(s[0:N])/param.norm
    jsx += 2.0 * drv * param.hx
    jsy = 2.0 * param.jy * param.jmat.dot(s[N:2*N])/param.norm
    jsy += 2.0 * drv * param.hy
    jsz = 2.0 * param.jz * param.jmat.dot(s[2*N:3*N])/param.norm
    jsz += 2.0 * drv * param.hz
    dsxdt = s[N:2*N] * jsz - s[2*N:3*N] * jsy
    dsydt = s[2*N:3*N] * jsx - s[0:N] * jsz
    dszdt = s[0:N] * jsy - s[N:2*N] * jsx
    return np.concatenate((dsxdt, dsydt, dszdt))

def jac_1storder(s, t, param):
    """
    Jacobian of the general case. First order.
    This is given by 9 NXN submatrices:
    J00=J11=J22=0
    Although Jacobian is NOT antisymmetric in general! See below
    J01 = +J_z diag(J|s^x>) + h(t) h_z - J_y (J#|s^z>)
    J10 = -J_z diag(J|s^x>) - h(t) h_z + J_x (J#|s^z>)
    J02 = -J_y diag(J|s^y>) - h(t) h_y + J_z (J#|s^y>)
    J20 = +J_y diag(J|s^y>) + h(t) h_y - J_x (J#|s^y>)
    J12 = +J_x diag(J|s^x>) + h(t) h_x - J_z (J#|s^x>)
    J21 = -J_x diag(J|s^x>) - h(t) h_x + J_y (J#|s^x>)
    Here, '#' (hash operator) means multiply each row of a matrix by the
    corresponding vector element. This is implemented by numpy.multiply()
    """
    N = param.latsize
    #s[0:N] = sx , s[N:2*N] = sy, s[2*N:3*N] = sz
    full_jacobian = np.zeros(shape=(3*N, 3*N))
    drivemat = 2.0 * drive(t, param) * np.eye(N)
    diag_jsx = np.diagflat((param.jmat.dot(s[0:N])))/param.norm
    diag_jsy = np.diagflat((param.jmat.dot(s[N:2*N])))/param.norm
    #diag_jsz = np.diagflat((param.jmat.dot(s[2*N:3*N])))/param.norm
    hash_jsx = (np.multiply(param.jmat.T, s[0:N]).T)/param.norm
    hash_jsy = (np.multiply(param.jmat.T, s[N:2*N]).T)/param.norm
    hash_jsz = (np.multiply(param.jmat.T, s[2*N:3*N]).T)/param.norm
    full_jacobian[0:N, N:2*N] = param.jz * diag_jsx + drivemat * param.hz\
      -param.jy * hash_jsz
    full_jacobian[N:2*N, 0:N] = -param.jz * diag_jsx - \
      drivemat * param.hz + param.jx * hash_jsz
    full_jacobian[0:N, 2*N:3*N] = -param.jy * diag_jsy - drivemat * \
      param.hy + param.jz * hash_jsy
    full_jacobian[2*N:3*N, 0:N] = param.jy * diag_jsy + drivemat * \
      param.hy - param.jx * hash_jsy
    full_jacobian[N:2*N, 2*N:3*N] = param.jx * diag_jsx + drivemat * \
      param.hx - param.jz * hash_jsx
    full_jacobian[2*N:3*N, N:2*N] = -param.jx * diag_jsx - drivemat * \
      param.hx + param.jy * hash_jsx
    return full_jacobian

def func_2ndorder(s, t, param):
    """
    The RHS of general case, second order correction, per Lorenzo
    "J" is the J_{ij} hopping matrix
    -\partial_t |s^x> = -first order + 2 (J^y  Jg^{yz} - J^z  Jg^{zy})
                                                                  /norm,
    -\partial_t |s^y> = -first order + 2 (-J^z  Jg^{zx} + J^x  Jg^{xz})
                                                                  /norm,
    -\partial_t |s^z> = -first order + 2 (-J^x  Jg^{xy} + J^y  Jg^{yx})
                                                                  /norm.
    """
    N = param.latsize
    #svec  is the tensor s^l_\mu
    #G = s[3*N:].reshape(3,3,N,N) is the tensor g^{ab}_{\mu\nu}.
    sview = s.view()
    stensor = sview[0:3*N].reshape(3, N)
    gtensor = sview[3*N:].reshape(3, 3, N, N)
    gtensor[:,:,range(N),range(N)] = 0.0 #Set the diagonals of g_munu to 0

    htensor = np.zeros_like(stensor)
    htensor[0].fill(param.hvec[0])
    htensor[1].fill(param.hvec[1])
    htensor[2].fill(param.hvec[2])
    Gtensor = np.einsum("mg,abgn->abmn", param.jmat, gtensor)/param.norm
    Mtensor = np.einsum("am,b,mn->abmn", stensor, param.jvec, \
      param.jmat)/param.norm
    hvec_dressed = htensor + np.einsum("llgm->lm", Mtensor)
    dtensor = gtensor + np.einsum("am,bn", stensor, stensor)
    dsdt_1 = func_1storder(sview[0:3*N], t, param).reshape(3, N)
    dsdt = dsdt_1 - \
      2.0 * np.einsum("bcmm,b,abc->am", Gtensor, param.jvec, eijk)

    dgdt = -np.einsum("lbmn,abl->abmn", Mtensor, eijk) + \
      np.einsum("lanm,abl->abmn", Mtensor, eijk)

    dgdt -= np.einsum("lm,kbmn,lka->abmn", hvec_dressed, gtensor, eijk) -\
      np.einsum("llnm,kbmn,lka->abmn", Mtensor, gtensor, eijk) +\
        np.einsum("ln,akmn,lkb->abmn", hvec_dressed, gtensor, eijk) -\
          np.einsum("llmn,akmn,lkb->abmn", Mtensor, gtensor, eijk)

    dgdt -= np.einsum("l,km,lbmn,lka->abmn", \
      param.jvec, stensor, Gtensor, eijk) + \
        np.einsum("l,kn,lanm,lkb->abmn", param.jvec, stensor, \
          Gtensor, eijk)

    dgdt += np.einsum("almn,lkmn,lkb->abmn", Mtensor, dtensor, eijk)\
      + np.einsum("blnm,lknm,lka->abmn", Mtensor, dtensor, eijk)

    #Flatten it before returning
    return np.concatenate((dsdt.flatten(), 2.0 * dgdt.flatten()))

def jac_2ndorder(s, t, param):
    """
    Jacobian of the general case. Second order.
    """
    N = param.latsize
    fullsize_2ndorder = 3 * N + 9 * N**2
    #svec  is the tensor s^l_\mu
    #G = s[3*N:].reshape(3,3,N,N) is the tensor g^{ab}_{\mu\nu}.
    sview = s.view()
    stensor = sview[0:3*N].reshape(3, N)
    gtensor = sview[3*N:].reshape(3, 3, N, N)
    htensor = np.zeros_like(stensor)
    htensor[0].fill(param.hvec[0])
    htensor[1].fill(param.hvec[1])
    htensor[2].fill(param.hvec[2])
    jjtensor = np.einsum("a,mn->amn", param.jvec, param.jmat)
    sstensor = np.einsum("km,ln->klmn",stensor,stensor)
    Mtensor = np.einsum("am,b,mn->abmn", stensor, param.jvec, \
          param.jmat)/param.norm
    hvec_dressed = htensor + np.einsum("llgm->lm", Mtensor)


    full_jacobian = np.zeros(shape=(fullsize_2ndorder, fullsize_2ndorder))

    #J00 subblock
    full_jacobian[0:3*N, 0:3*N] = jac_1storder(s, t, param)

    #J01 subblock. Precalculated
    full_jacobian[0:3*N, 3*N:] = param.dsdotdg

    #J10 subblock
    full_jacobian[3*N:, 0:3*N] =  -(np.einsum("pml,kbmn,pka->abpmnl", \
      jjtensor,gtensor, eijk) +  np.einsum("pnl,akmn,pkb->abpmnl", \
        jjtensor, gtensor, eijk)).reshape(9*N*N,3*N)
    full_jacobian[3*N:, 0:3*N] -= (np.einsum("qmg,ml,bqng,qpa->abpmnl",\
      jjtensor, param.deltamn,gtensor, eijk) + \
        np.einsum("qng,nl,aqmg,qpb->abpmnl",jjtensor, param.deltamn, \
          gtensor, eijk) ).reshape(9*N*N,3*N)
    full_jacobian[3*N:, 0:3*N] += (np.einsum("qmn,ml,bqnn,qpa->abpmnl",\
      jjtensor, param.deltamn,gtensor, eijk) + \
        np.einsum("qnm,nl,aqmm,qpb->abpmnl", jjtensor,param.deltamn, \
          gtensor, eijk)).reshape(9*N*N,3*N)
    full_jacobian[3*N:, 0:3*N] += (np.einsum("qmn,ml,pa,qkmn,qkb->abpmnl",\
      jjtensor,param.deltamn,deltaij,gtensor+sstensor,eijk) + \
        np.einsum("qmn,nl,pb,kqmn,qka->abpmnl", jjtensor,param.deltamn, \
          deltaij,gtensor+sstensor,eijk)).reshape(9*N*N,3*N)
    full_jacobian[3*N:, 0:3*N] += (np.einsum("pmn,ml,akmn,pkb->abpmnl",\
      jjtensor,param.deltamn, sstensor, eijk) + \
        np.einsum("pmn,nl,bknm,pka->abpmnl", jjtensor,param.deltamn, \
          sstensor, eijk) + np.einsum("kmn,nl,akmm,kpb->abpmnl",\
            jjtensor,param.deltamn, sstensor, eijk) + \
              np.einsum("kmn,ml,bknn,kpa->abpmnl", jjtensor,param.deltamn, \
                sstensor, eijk)).reshape(9*N*N,3*N)
    full_jacobian[3*N:, 0:3*N] = 2.0 * \
      (full_jacobian[3*N:, 0:3*N]/param.norm)
    full_jacobian[3*N:, 0:3*N] += param.dsdotdg.T

    #J11 subblock:
    full_jacobian[3*N:, 3*N:] = -(np.einsum("qm,mlnhbpqra->abrpmnlh",\
       hvec_dressed, param.delta_eps_tensor)).reshape(9*N*N,9*N*N)
    full_jacobian[3*N:, 3*N:] += (np.einsum("qqmn,mlnhbpqra->abrpmnlh", \
          Mtensor, param.delta_eps_tensor)).reshape(9*N*N,9*N*N)
    full_jacobian[3*N:, 3*N:] -= (np.einsum("qn,mlnharqpb->abrpmnlh",\
      hvec_dressed, param.delta_eps_tensor)).reshape(9*N*N,9*N*N)
    full_jacobian[3*N:, 3*N:] += (np.einsum("qqnm,mlnharqpb->abrpmnlh",\
          Mtensor, param.delta_eps_tensor)).reshape(9*N*N,9*N*N)

    excl_tensor  = -np.einsum("qmh,km,nl,br,pka->abrpmnlh",\
          jjtensor,stensor, param.deltamn, deltaij,eijk)
    excl_tensor += -np.einsum("qnh,kn,ml,ar,pkb->abrpmnlh",\
          jjtensor,stensor, param.deltamn, deltaij,eijk)
    excl_tensor += -np.einsum("qml,km,nh,bp,rka->abrpmnlh",\
          jjtensor,stensor, param.deltamn, deltaij,eijk)
    excl_tensor += -np.einsum("qnl,kn,mh,ap,rkb->abrpmnlh",\
          jjtensor,stensor, param.deltamn, deltaij,eijk)
    #Set the \eta=\mu,\nu components of excl_tensor to 0
    excl_tensor[:,:,:,:,range(N),:,:,range(N)] = 0.0
    excl_tensor[:,:,:,:,:,range(N),:,range(N)] = 0.0
    full_jacobian[3*N:, 3*N:] += excl_tensor.reshape(9*N*N,9*N*N)

    full_jacobian[3*N:, 3*N:] += (np.einsum("rmn,am,ml,nh,rpb->abrpmnlh",\
      jjtensor,stensor,param.deltamn,param.deltamn,eijk) + \
        np.einsum("rmn,bn,mh,nl,rpa->abrpmnlh",\
          jjtensor,stensor,param.deltamn,param.deltamn,eijk)).reshape(9*N*N,9*N*N)
    full_jacobian[3*N:, 3*N:] -= (np.einsum("pmn,am,mh,nl,prb->abrpmnlh",\
      jjtensor,stensor,param.deltamn,param.deltamn,eijk) + \
        np.einsum("pmn,bn,ml,nh,pra->abrpmnlh",\
          jjtensor,stensor,param.deltamn,param.deltamn,eijk)).reshape(9*N*N,9*N*N)
    full_jacobian[3*N:, 3*N:] = 2.0 * (full_jacobian[3*N:, 3*N:]/param.norm)

    return full_jacobian

class ParamData:
    """Class to store parameters, precalculated objects,
        filenames, objects like Kac norm
        time-independent part of Jacobian. Set s_order to
        true if doing second order dtwa'
    """
    def __init__(self,pbc=False ,nonorm=True, latsize=101, beta=1.0, \
                          h0=1.0, omega=0.0, hx=0.0, hy=0.0, hz=0.0,\
                            jx=0.0, jy=0.0, jz=1.0):
        #Default Output file names. Each file dumps a different observable
        self.output_magx = "sx_outfile.txt"
        self.output_magy = "sy_outfile.txt"
        self.output_magz = "sz_outfile.txt"

        self.output_sxvar = "sxvar_outfile.txt"
        self.output_syvar = "syvar_outfile.txt"
        self.output_szvar = "szvar_outfile.txt"

        self.output_sxyvar = "sxyvar_outfile.txt"
        self.output_sxzvar = "sxzvar_outfile.txt"
        self.output_syzvar = "syzvar_outfile.txt"

        #Whether to normalize with Kac norm or not
        self.nonorm = nonorm

        self.latsize = latsize

        self.beta = beta #Power law index for long range interactions
        self.h0 = h0 # Drive amplitude
        self.omega = omega #Drive frequency
        self.hx = hx #x transverse field
        self.hy = hy #y transverse field
        self.hz = hz #z transverse field
        self.jx = jx #x hopping
        self.jy = jy #y hopping
        self.jz = jz #z hopping

        self.jvec = np.array([jx, jy, jz])
        self.hvec = np.array([hx, hy, hz])
        N = self.latsize
        self.fullsize_2ndorder = 3 * N + 9 * N**2
        self.deltamn = np.eye(N)
        # These are the lattice  sites for two point density matrix calc.
        self.tpnt_sites = (np.floor(N/2), np.floor(N/2)+2)
        if(pbc):
            self.periodic_boundary_conditions = True
            self.open_boundary_conditions = False
            #This is the dense Jmn hopping matrix with power law decay for
            #periodic or open boundary conditions.
            J = dia_matrix((N, N))
            mid_diag = np.floor(N/2).astype(int)
            for i in xrange(1,mid_diag+1):
                elem = pow(i, -self.beta)
                J.setdiag(elem, k=i)
                J.setdiag(elem, k=-i)
            for i in xrange(mid_diag+1, N):
                elem = pow(N-i, -self.beta)
                J.setdiag(elem, k=i)
                J.setdiag(elem, k=-i)
        else: #Open boundary conditions
            self.periodic_boundary_conditions = False
            self.open_boundary_conditions = True
            J = dia_matrix((N, N))
            for i in xrange(1,N):
                elem = pow(i, -self.beta)
                J.setdiag(elem, k=i)
                J.setdiag(elem, k=-i)

        self.jmat = J.toarray()

        #This is the optional Kac norm
        mid = np.floor(N/2).astype(int)
        if self.nonorm:
            self.norm = 1.0
        else:
            self.norm =\
                  2.0 * np.sum(1/(pow(\
                    np.arange(1, mid+1), self.beta).astype(float)))

class OutData:
    """Class to store output data"""
    def __init__(self, t, sx, sy, sz, sxx, syy, szz, sxy, sxz, syz,\
      params):
        self.t_output = t
        self.sx, self.sy, self.sz = sx, sy, sz
        self.sxvar, self.syvar, self.szvar = sxx, syy, szz
        self.sxyvar, self.sxzvar, self.syzvar = sxy, sxz, syz
        self.__dict__.update(params.__dict__)

    def normalize_data(self, w_totals, lsize):
        self.sx = self.sx/(w_totals * lsize)
        self.sy = self.sy/(w_totals * lsize)
        self.sz = self.sz/(w_totals * lsize)
        self.sxvar = (1/lsize) + (self.sxvar/(w_totals * lsize * lsize))
        self.sxvar = self.sxvar - (self.sx)**2
        self.syvar = (1/lsize) + (self.syvar/(w_totals * lsize * lsize))
        self.syvar = self.syvar - (self.sy)**2
        self.szvar = (1/lsize) + (self.szvar/(w_totals * lsize * lsize))
        self.szvar = self.szvar - (self.sz)**2
        self.sxyvar = (self.sxyvar/(w_totals * lsize * lsize))
        self.sxyvar = self.sxyvar - (self.sx * self.sy)
        self.sxzvar = (self.sxzvar/(w_totals * lsize * lsize))
        self.sxzvar = self.sxzvar - (self.sx * self.sz)
        self.syzvar = (self.syzvar/(w_totals * lsize * lsize))
        self.syzvar = self.syzvar - (self.sy * self.sz)

    def dump_data(self):
        np.savetxt(self.output_magx, \
          np.vstack((self.t_output, self.sx)).T, delimiter=' ')
        np.savetxt(self.output_magy, \
          np.vstack((self.t_output, self.sy)).T, delimiter=' ')
        np.savetxt(self.output_magz, \
          np.vstack((self.t_output, self.sz)).T, delimiter=' ')
        np.savetxt(self.output_sxvar, \
          np.vstack((self.t_output, self.sxvar)).T, delimiter=' ')
        np.savetxt(self.output_syvar, \
          np.vstack((self.t_output, self.syvar)).T, delimiter=' ')
        np.savetxt(self.output_szvar, \
          np.vstack((self.t_output, self.szvar)).T, delimiter=' ')
        np.savetxt(self.output_sxyvar, \
          np.vstack((self.t_output, self.sxyvar)).T, delimiter=' ')
        np.savetxt(self.output_sxzvar, \
          np.vstack((self.t_output, self.sxzvar)).T, delimiter=' ')
        np.savetxt(self.output_syzvar, \
          np.vstack((self.t_output, self.syzvar)).T, delimiter=' ')

class OutData_ij:
    """
    Class to store ij output data
    The gij is a numpy array of 3X3 gij matrices at multiple times
    """
    def __init__(self, t, sites, sxi, syi, szi, sxj, syj, szj, sy_iplusk,\
      syy_k, gij=None):
        self.times = t
        self.sites = sites
        self.sxi, self.syi, self.szi = sxi, syi, szi
        self.sxj, self.syj, self.szj = sxj, syj, szj
        #Output formatting dictionaries
        self.sitespinsdict = {"time": t,\
                "sxi": self.sxi.view(),\
                "syi": self.syi.view(),\
                "szi": self.szi.view(),\
                "sxj": self.sxj.view(),\
                "syj": self.syj.view(),\
                "szj": self.szj.view()}
        self.sy_iplusk = sy_iplusk
        self.syy_k = syy_k
        if gij is not None:
            self.gij = gij
            v = self.gij.view()
            self.sitecorrdict = {"time": t,\
                    "gxxij": v[:, 0, 0],\
                    "gxyij": v[:, 0, 1],\
                    "gxzij": v[:, 0, 2],\
                    "gyxij": v[:, 1, 0],\
                    "gyyij": v[:, 1, 1],\
                    "gyzij": v[:, 1, 2],\
                    "gzxij": v[:, 2, 0],\
                    "gzyij": v[:, 2, 1],\
                    "gzzij": v[:, 2, 2]}

    def normalize_data(self, w_totals):
        self.sitespinsdict['sxi'] = self.sitespinsdict['sxi']/(w_totals)
        self.sitespinsdict['syi'] = self.sitespinsdict['syi']/(w_totals)
        self.sitespinsdict['szi'] = self.sitespinsdict['szi']/(w_totals)
        self.sitespinsdict['sxj'] = self.sitespinsdict['sxj']/(w_totals)
        self.sitespinsdict['syj'] = self.sitespinsdict['syj']/(w_totals)
        self.sitespinsdict['szj'] = self.sitespinsdict['szj']/(w_totals)
        #Normalize the spatial correlations:
        self.sy_iplusk = self.sy_iplusk/(w_totals)
        self.syy_k = self.syy_k/(w_totals)
        self.syy_k -= np.array([self.sy_iplusk[i] *\
          self.sitespinsdict['syi'][i] for i in xrange(self.times.size)])

        if hasattr(self, 'sitecorrdict'):
            for key in self.sitecorrdict.iterkeys():
                if key is not "time":
                    self.sitecorrdict[key] =  self.sitecorrdict[key]/(w_totals)

    def dump_data(self):
        print("\n\n Tabular dump of site data:\n\n")
        print("\n Note that, in the first order case,")
        print("the 'gij' columns actually print sijs,")
        print("which are s(..)i * s(..)j\n\n")
        print("Sites chosen:", self.sites)
        print("\n")
        print(tabulate(self.sitespinsdict, headers="keys", floatfmt=".6f" ))
        print("Spatial correlations from site i = \n", self.sites[0])
        print(np.vstack((self.times,self.syy_k.T)).T)

        if hasattr(self, 'sitecorrdict'):
            print("    ")
            print(tabulate(self.sitecorrdict, headers="keys", floatfmt=".6f"))

class Dtwa_System:
    """
    This is the class that creates the DTWA system,
    has all MPI_Gather routines for aggregating the
    samples, and executes the dtwa methods (1st and 2nd order)
    Set s_order to true if doing second order dtwa
    Set jac to false if you don't want to evaluate the jacobian, since
    it may be too big in some cases and cause the routine to crash.
    """

    def __init__(self, params, mpicomm, n_t=2000, file_output=True,\
                            seed_offset=0,  s_order=False, jac=False,\
                              verbose=True, sitedata=False):
        """
        Input default values and amke precalculated objects
        Comm = MPI Communicator
        """
        self.jac = jac
        self.__dict__.update(params.__dict__)
        self.n_t = n_t
        self.file_output = file_output
        self.comm=mpicomm
        self.seed_offset = seed_offset
        self.s_order = s_order
        #Booleans for verbosity and for calculating site data
        self.verbose = verbose
        self.sitedata = sitedata
        N = params.latsize

        #Only computes these if you want 2nd order
        if self.s_order and self.jac:
            #Below are the constant subblocks of the 2nd order Jacobian
            #The 00 subblock is the first order Jacobian in func below
            #The entire 01 subblock, fully time independent (ds_dot/dg):
            self.dsdotdg = -np.einsum("p,mh,ml,apr->arpmlh",\
              self.jvec, self.jmat, self.deltamn, eijk)
            self.dsdotdg += np.einsum("r,ml,mh,arp->arpmlh", \
              self.jvec,self.jmat, self.deltamn, eijk)
            self.dsdotdg = 2.0 * (self.dsdotdg/self.norm)
            self.dsdotdg = self.dsdotdg.reshape(3*N, 9*N**2)
            self.delta_eps_tensor = np.einsum("ml,nh,ar,qpb->mlnharqpb",\
              self.deltamn,self.deltamn,deltaij,eijk)
            self.delta_eps_tensor += np.einsum("mh,nl,ap,qrb->mhnlapqrb",\
              self.deltamn,self.deltamn,deltaij,eijk)
            #The time independent part of the 10 subblock (dg_dot/ds):
            #is the SAME as ds_dot/dg

    def sum_reduce_all_data(self, datalist_loc,t, mpcomm):
        """
        Does the parallel sum reduction of all data
        """
        #Do local sums
        sx_locsum = np.sum(data.sx for data in datalist_loc)
        sy_locsum = np.sum(data.sy for data in datalist_loc)
        sz_locsum = np.sum(data.sz for data in datalist_loc)
        sxvar_locsum = np.sum(data.sxvar for data in datalist_loc)
        syvar_locsum = np.sum(data.syvar for data in datalist_loc)
        szvar_locsum = np.sum(data.szvar for data in datalist_loc)
        sxyvar_locsum = np.sum(data.sxyvar for data in datalist_loc)
        sxzvar_locsum = np.sum(data.sxzvar for data in datalist_loc)
        syzvar_locsum = np.sum(data.syzvar for data in datalist_loc)

        #Only root processor will actually get the data
        sx_totals = np.zeros_like(sx_locsum) if mpcomm.rank == root\
          else None
        sy_totals = np.zeros_like(sy_locsum) if mpcomm.rank == root\
          else None
        sz_totals = np.zeros_like(sz_locsum) if mpcomm.rank == root\
          else None
        sxvar_totals = np.zeros_like(sxvar_locsum) if mpcomm.rank == root\
          else None
        syvar_totals = np.zeros_like(syvar_locsum) if mpcomm.rank == root\
          else None
        szvar_totals = np.zeros_like(szvar_locsum) if mpcomm.rank == root\
          else None
        sxyvar_totals = np.zeros_like(sxyvar_locsum) if mpcomm.rank == root\
          else None
        sxzvar_totals = np.zeros_like(sxzvar_locsum) if mpcomm.rank == root\
          else None
        syzvar_totals = np.zeros_like(syzvar_locsum) if mpcomm.rank == root\
          else None

        #To prevent conflicts with other comms
        duplicate_comm = Intracomm(mpcomm)
        sx_totals = duplicate_comm.reduce(sx_locsum, root=root)
        sy_totals = duplicate_comm.reduce(sy_locsum, root=root)
        sz_totals = duplicate_comm.reduce(sz_locsum, root=root)
        sxvar_totals = duplicate_comm.reduce(sxvar_locsum, root=root)
        syvar_totals = duplicate_comm.reduce(syvar_locsum, root=root)
        szvar_totals = duplicate_comm.reduce(szvar_locsum, root=root)
        sxyvar_totals = duplicate_comm.reduce(sxyvar_locsum, root=root)
        sxzvar_totals = duplicate_comm.reduce(sxzvar_locsum, root=root)
        syzvar_totals = duplicate_comm.reduce(syzvar_locsum, root=root)

        if mpcomm.rank == root:
            return OutData(t, sx_totals, sy_totals, sz_totals, sxvar_totals, \
                syvar_totals, szvar_totals, sxyvar_totals, sxzvar_totals,\
                  syzvar_totals, self)
        else:
            return None

    def sum_reduce_site_data(self, datalist_loc, t, sites, mpcomm):
        """
        Does the parallel sum reduction of site data
        """
        sxi_locsum = np.sum(data.sxi for data in datalist_loc)
        syi_locsum = np.sum(data.syi for data in datalist_loc)
        szi_locsum = np.sum(data.szi for data in datalist_loc)
        sxj_locsum = np.sum(data.sxj for data in datalist_loc)
        syj_locsum = np.sum(data.syj for data in datalist_loc)
        szj_locsum = np.sum(data.szj for data in datalist_loc)
        sy_iplusk_locsum = np.sum(data.sy_iplusk for data in datalist_loc)
        syy_k_locsum = np.sum(data.syy_k for data in datalist_loc)

        try: #This is to take care of the case when gij = None
            gijs_locsum = np.sum(data.gij for data in datalist_loc)
        except AttributeError:
            gijs_locsum = None

        sxi_totals = np.zeros_like(sxi_locsum) if mpcomm.rank == root\
          else None
        syi_totals = np.zeros_like(syi_locsum) if mpcomm.rank == root\
          else None
        szi_totals = np.zeros_like(szi_locsum) if mpcomm.rank == root\
          else None
        sxj_totals = np.zeros_like(sxj_locsum) if mpcomm.rank == root\
          else None
        syj_totals = np.zeros_like(syj_locsum) if mpcomm.rank == root \
          else None
        szj_totals = np.zeros_like(szj_locsum) if mpcomm.rank == root \
          else None
        sy_iplusk_totals = np.zeros_like(sy_iplusk_locsum) \
          if mpcomm.rank == root else None
        syy_k_totals = np.zeros_like(syy_k_locsum) \
          if mpcomm.rank == root else None

        gijs_totals = np.zeros_like(gijs_locsum) if mpcomm.rank == root \
          else None

        #To prevent conflicts with other comms
        duplicate_comm = Intracomm(mpcomm)
        #Broadcast these reductions to root
        sxi_totals = duplicate_comm.reduce(sxi_locsum, root=root)
        syi_totals = duplicate_comm.reduce(syi_locsum, root=root)
        szi_totals = duplicate_comm.reduce(szi_locsum, root=root)
        sxj_totals = duplicate_comm.reduce(sxj_locsum, root=root)
        syj_totals = duplicate_comm.reduce(syj_locsum, root=root)
        szj_totals = duplicate_comm.reduce(szj_locsum, root=root)
        sy_iplusk_totals = duplicate_comm.reduce(sy_iplusk_locsum,root=root)
        syy_k_totals = duplicate_comm.reduce(syy_k_locsum, root=root)

        if gijs_locsum is not None:
            gijs_totals = duplicate_comm.reduce(gijs_locsum, root=root)
        else:
            gijs_totals = None

        if mpcomm.rank == root:
            return OutData_ij(t, sites, sxi_totals, syi_totals, \
              szi_totals, sxj_totals, syj_totals, szj_totals, \
                sy_iplusk_totals, syy_k_totals, gijs_totals)
        else:
            return None

    def dtwa_ising_longrange_1storder(self, time_info):
        comm = self.comm
        N = self.latsize
        (t_init, n_cycles, n_steps) = time_info
        rank = comm.rank

        if rank == root and self.verbose:
            pprint("# Run parameters:")
            pprint(vars(self), depth=2)
        if rank == root and not self.verbose:
            pprint("# Starting run ...")
        if self.omega == 0:
            t_final = t_init + n_cycles
        else:
            t_final = t_init + (n_cycles * (2.0* np.pi/self.omega))
        dt = (t_final-t_init)/(n_steps-1.0)
        t_output = np.arange(t_init, t_final, dt)
        #Let each process get its chunk of n_t by round robin
        nt_loc = 0
        iterator = rank
        while iterator < self.n_t:
            nt_loc += 1
            iterator += comm.size
        #Scatter unique seeds for generating unique random number arrays :
        #each processor gets its own nt_loc seeds, and allocates nt_loc
        #initial conditions. Each i.c. is a 2N sized array
        #now, each process sends its value of nt_loc to root
        all_ntlocs = comm.gather(nt_loc, root=root)
        #Let the root process initialize nt unique integers for random seeds
        if rank == root:
            all_seeds = np.arange(self.n_t, dtype=np.int64)+1
            all_ntlocs = np.array(all_ntlocs)
            all_displacements = np.roll(np.cumsum(all_ntlocs), root+1)
            all_displacements[root] = 0 # First displacement
        else:
            all_seeds = None
            all_displacements = None
        local_seeds = np.zeros(nt_loc, dtype=np.int64)
        #Root scatters nt_loc sized seed data to that particular process
        comm.Scatterv([all_seeds, all_ntlocs, all_displacements,\
          MPI.DOUBLE], local_seeds)

        list_of_local_data = []

        if self.sitedata:
            list_of_local_ijdata = []

        for runcount in xrange(0, nt_loc, 1):
            random.seed(local_seeds[runcount] + self.seed_offset)
            #According to Schachenmayer, the wigner function of the quantum
            #state generates the below initial conditions classically
            sx_init = np.ones(N)
            sy_init = 2.0 * np.random.randint(0,2, size=N) - 1.0
            sz_init = 2.0 * np.random.randint(0,2, size=N) - 1.0
            #Set initial conditions for the dynamics locally to vector
            #s_init and store it as [s^x,s^x,s^x, .... s^y,s^y,s^y ...,
            #s^z,s^z,s^z, ...]
            s_init = np.concatenate((sx_init, sy_init, sz_init))
            if self.verbose:
                if self.jac:
                    s, info = odeint(func_1storder, s_init, t_output,\
                      args=(self,), Dfun=jac_1storder, full_output=True)
                else:
                    s, info = odeint(func_1storder, s_init, t_output,\
                      args=(self,), Dfun=None, full_output=True)
            else:
                if self.jac:
                    s = odeint(func_1storder, s_init, t_output, args=(self,),\
                      Dfun=jac_1storder)
                else:
                    s = odeint(func_1storder, s_init, t_output, args=(self,),\
                      Dfun=None)
            #Compute expectations <sx> and \sum_{ij}<sx_i sx_j> -<sx>^2 with
            #wigner func at t_output values LOCALLY for each initcond and
            #store them
            sx_expectations = np.sum(s[:, 0:N], axis=1)
            sy_expectations = np.sum(s[:, N:2*N], axis=1)
            sz_expectations = np.sum(s[:, 2*N:3*N], axis=1)

            if self.sitedata:
                (i, j) = self.tpnt_sites
                sxi, syi, szi = s[:, i], s[:, i+N], s[:, i+2*N]
                sxj, syj, szj = s[:, j], s[:, j+N], s[:, j+2*N]
                sxxij, syyij, szzij = sxi * sxj, syi * syj, szi * szj
                sxyij, sxzij, syzij = sxi * syj, sxi * szj, syi * szj
                syxij, szxij, szyij = syi * sxj, szi * sxj, szi * syj
                gij = np.array([sxxij, sxyij, sxzij, syxij, syyij, syzij,\
                  szxij, szyij, szzij]).T.reshape(t_output.size,3,3)
                sxi, syi, szi = sxi , syi , \
                  szi
                sxj, syj, szj = sxj , syj , \
                  szj
                #Calculate Spatial Correlations
                sy_iplusk = s[:, N:2*N][:,i:] #This is a matrix
                syy_k = np.array([sy_iplusk[t] * syi[t] \
                  for t in xrange(t_output.size)])# This is also a matrix

                localdataij = OutData_ij(t_output, self.tpnt_sites, \
                                                          sxi, syi, szi,\
                                                            sxj, syj, szj,\
                                                                sy_iplusk,\
                                                                  syy_k,\
                                                                        gij)
                list_of_local_ijdata.append(localdataij)

            #Quantum spin variance maps to the classical expression
            # (1/N) + (1/N^2)\sum_{i\neq j} S^x_i S^x_j - <S^x>^2 and
            # (1/N) + (1/N^2)\sum_{i\neq j} S^y_i S^z_j
            # since the i=j terms quantum average to unity
            sx_var =   (np.sum(s[:, 0:N], axis=1)**2 \
              - np.sum(s[:, 0:N]**2, axis=1))
            sy_var =   (np.sum(s[:, N:2*N], axis=1)**2 \
            - np.sum(s[:, N:2*N]**2, axis=1))
            sz_var =   (np.sum(s[:, 2*N:3*N], axis=1)**2 \
            - np.sum(s[:, 2*N:3*N]**2, axis=1))

            sxy_var =   np.sum([fftconvolve(s[m, 0:N], \
              s[m, N:2*N]) for m in xrange(t_output.size)], axis=1)
            sxz_var =   np.sum([fftconvolve(s[m, 0:N], \
              s[m, 2*N:3*N]) for m in xrange(t_output.size)], axis=1)
            syz_var =   np.sum([fftconvolve(s[m, N:2*N], \
              s[m, 2*N:3*N]) for m in xrange(t_output.size)], axis=1)

            localdata = OutData(t_output, sx_expectations, sy_expectations,\
              sz_expectations, sx_var, sy_var, sz_var, sxy_var, sxz_var, \
                syz_var, self)
            list_of_local_data.append(localdata)
        #After loop above  sum reduce (don't forget to average) all locally
        #calculated expectations at each time to root
        outdat = \
          self.sum_reduce_all_data(list_of_local_data, t_output, comm)
        if self.sitedata:
            sij = self.sum_reduce_site_data(list_of_local_ijdata,\
                                            t_output, self.tpnt_sites, comm)
            if rank == root:
                sij.normalize_data(self.n_t)
                if self.file_output:
                    sij.dump_data()

        if rank == root:
                #Dump to file
            outdat.normalize_data(self.n_t, N)
            if self.file_output:
                outdat.dump_data()
            if self.verbose:
                print("  ")
                print("Integration output info:")
                pprint(info)
                print("""# Cumulative number of Jacobian evaluations
                            by root:""", \
                  np.sum(info['nje']))
            print('# Done!')
            return outdat
        else:
            return None

    def dtwa_ising_longrange_2ndorder(self, time_info, sampling):
        old_settings = np.seterr(all='ignore') #Prevent overflow warnings
        comm=self.comm
        N = self.latsize
        (t_init, n_cycles, n_steps) = time_info
        rank = comm.rank
        if rank == root and self.verbose:
            pprint("# Run parameters:")
            #Copy params to another object, then delete
            #the output that you don't want printed
            out = copy.copy(self)
            out.dsdotdg = 0.0
            out.delta_eps_tensor = 0.0
            out.jmat = 0.0
            out.deltamn = 0.0
            pprint(vars(out), depth=2)
        if rank == root and not self.verbose:
            pprint("# Starting run ...")
        if self.omega == 0:
            t_final = t_init + n_cycles
        else:
            t_final = t_init + (n_cycles * (2.0* np.pi/self.omega))
        dt = (t_final-t_init)/(n_steps-1.0)
        t_output = np.arange(t_init, t_final, dt)
        #Let each process get its chunk of n_t by round robin
        nt_loc = 0
        iterator = rank
        while iterator < self.n_t:
            nt_loc += 1
            iterator += comm.size
        #Scatter unique seeds for generating unique random number arrays :
        #each processor gets its own nt_loc seeds, and allocates nt_loc
        #initial conditions. Each i.c. is a 2N sized array
        #now, each process sends its value of nt_loc to root
        all_ntlocs = comm.gather(nt_loc, root=root)
        #Let the root process initialize nt unique integers for random seeds
        if rank == root:
            all_seeds = np.arange(self.n_t, dtype=np.int64)+1
            all_ntlocs = np.array(all_ntlocs)
            all_displacements = np.roll(np.cumsum(all_ntlocs), root+1)
            all_displacements[root] = 0 # First displacement
        else:
            all_seeds = None
            all_displacements = None
        local_seeds = np.zeros(nt_loc, dtype=np.int64)
        #Root scatters nt_loc sized seed data to that particular process
        comm.Scatterv([all_seeds, all_ntlocs, all_displacements,\
          MPI.DOUBLE],local_seeds)

        list_of_local_data = []

        if self.verbose:
            list_of_dhwdt_abs2 = []

        if self.sitedata:
            list_of_local_ijdata = []

        for runcount in xrange(0, nt_loc, 1):
            random.seed(local_seeds[runcount] + self.seed_offset)
            sx_init = np.ones(N)
            if sampling == "spr":
                #According to Schachenmayer, the wigner function of the quantum
                #state generates the below initial conditions classically
                sy_init = 2.0 * np.random.randint(0,2, size=N) - 1.0
                sz_init = 2.0 * np.random.randint(0,2, size=N) - 1.0
                #Set initial conditions for the dynamics locally to vector
                #s_init and store it as [s^x,s^x,s^x, .... s^y,s^y,s^y ...,
                #s^z,s^z,s^z, ...]
                s_init_spins = np.concatenate((sx_init, sy_init, sz_init))
            elif sampling == "1-0":
                spin_choices = np.array([(1, 1,0),(1, 0,1),(1, -1,0),(1, 0,-1)])
                spins = np.array([random.choice(spin_choices) for i in xrange(N)])
                s_init_spins = spins.T.flatten()
            elif sampling == "all":
                spin_choices_spr = np.array([(1, 1,1),(1, 1,-1),(1, -1,1),(1, -1,-1)])
                spin_choices_10 = np.array([(1, 1,0),(1, 0,1),(1, -1,0),(1, 0,-1)])
                spin_choices = np.concatenate((spin_choices_10, spin_choices_spr))
                spins = np.array([random.choice(spin_choices) for i in xrange(N)])
                s_init_spins = spins.T.flatten()
            else:
                pass

            # Set initial correlations to 0.
            s_init_corrs = np.zeros(9*N*N)
            #Redirect unwanted stdout warning messages to /dev/null
            with stdout_redirected():
                if self.verbose:
                    if self.jac:
                        s, info = odeint(func_2ndorder, \
                          np.concatenate((s_init_spins, s_init_corrs)), t_output, \
                            args=(self,), Dfun=jac_2ndorder, full_output=True)
                    else:
                        s, info = odeint(func_2ndorder, \
                          np.concatenate((s_init_spins, s_init_corrs)),t_output, \
                            args=(self,), Dfun=None, full_output=True)
                else:
                    if self.jac:
                        s = odeint(func_2ndorder, \
                          np.concatenate((s_init_spins, s_init_corrs)), \
                            t_output, args=(self,), Dfun=jac_2ndorder)
                    else:
                        s = odeint(func_2ndorder, \
                          np.concatenate((s_init_spins, s_init_corrs)), t_output, \
                            args=(self,), Dfun=None)

            #Computes |dH/dt|^2 for a particular alphavec & weighes it
            #If the rms over alphavec of these are 0, then each H is const
            if self.verbose:
                hws = weyl_hamilt(s,t_output, self)
                dhwdt = np.array([t_deriv(hw, t_output) for hw in hws])
                dhwdt_abs2 = np.square(dhwdt)
                list_of_dhwdt_abs2.extend(dhwdt_abs2)

            s = np.array(s, dtype="float128")#Widen memory to reduce overflows
            #Compute expectations <sx> and \sum_{ij}<sx_i sx_j> -<sx>^2 with
            #wigner func at t_output values LOCALLY for each initcond and
            #store them
            sx_expectations = np.sum(s[:, 0:N], axis=1)
            sy_expectations = np.sum(s[:, N:2*N], axis=1)
            sz_expectations = np.sum(s[:, 2*N:3*N], axis=1)

            if self.sitedata:
                (i, j) = self.tpnt_sites
                sxi, syi, szi = s[:, i], s[:, i+N], s[:, i+2*N]
                sxi, syi, szi = sxi , syi ,\
                  szi
                sxj, syj, szj = s[:, j], s[:, j+N], s[:, j+2*N]
                sxj, syj, szj = sxj , syj ,\
                  szj
                sview = s.view()
                gij = sview[:,3*N:].reshape(\
                  t_output.size,3, 3, N, N)[:, :, :, i, j]
                #Calculate Spatial Correlations
                sy_iplusk = s[:, N:2*N][:,i:] #This is a matrix
                syy_k = np.array([sy_iplusk[t] * syi[t] \
                  for t in xrange(t_output.size)])# This is also a matrix
                localdataij = OutData_ij(t_output, self.tpnt_sites, \
                                                          sxi, syi, szi,\
                                                            sxj, syj, szj,\
                                                                sy_iplusk,\
                                                                  syy_k,\
                                                                        gij)
                list_of_local_ijdata.append(localdataij)

            #svec  is the tensor s^l_\mu
            #G = s[3*N:].reshape(3,3,N,N) is the tensor g^{ab}_{\mu\nu}.
            s = np.array(s, dtype="float128")#Enlarge in mem
            sview = s.view()
            gt = sview[:, 3*N:].reshape(s.shape[0], 3, 3, N, N)
            gt[:,:,:,range(N),range(N)] = 0.0 #Set diags to 0
            #Quantum spin variance
            sx_var = np.sum(gt[:,0,0,:,:], axis=(-1,-2))
            sx_var += (np.sum(s[:, 0:N], axis=1)**2 \
              - np.sum(s[:, 0:N]**2, axis=1))

            sy_var = np.sum(gt[:,1,1,:,:], axis=(-1,-2))
            sy_var += (np.sum(s[:, N:2*N], axis=1)**2 \
            - np.sum(s[:, N:2*N]**2, axis=1))

            sz_var = np.sum(gt[:,2,2,:,:], axis=(-1,-2))
            sz_var += (np.sum(s[:, 2*N:3*N], axis=1)**2 \
            - np.sum(s[:, 2*N:3*N]**2, axis=1))

            sxy_var = np.sum(gt[:,0,1,:,:], axis=(-1,-2))
            sxy_var += np.sum([fftconvolve(s[m, 0:N], s[m, N:2*N]) \
            for m in xrange(t_output.size)], axis=1)
            #Remove the diagonal parts
            sxy_var -= np.sum(s[:, 0:N] *  s[:, N:2*N], axis=1)

            sxz_var = np.sum(gt[:,0,2,:,:], axis=(-1,-2))
            sxz_var += np.sum([fftconvolve(s[m, 0:N], s[m, 2*N:3*N]) \
              for m in xrange(t_output.size)], axis=1)
            #Remove the diagonal parts
            sxz_var -= np.sum(s[:, 0:N] *  s[:, 2*N:3*N], axis=1)

            syz_var = np.sum(gt[:,1,2,:,:], axis=(-1,-2))
            syz_var += np.sum([fftconvolve(s[m, N:2*N], s[m, 2*N:3*N]) \
              for m in xrange(t_output.size)], axis=1)
            #Remove the diagonal parts
            syz_var -= np.sum(s[:, N:2*N] *  s[:, 2*N:3*N], axis=1)

            localdata = OutData(t_output, sx_expectations, sy_expectations,\
              sz_expectations, sx_var, sy_var, sz_var, sxy_var, sxz_var, \
                syz_var, self)
            list_of_local_data.append(localdata)
        #After loop above  sum reduce (don't forget to average) all locally
        #calculated expectations at each time to root
        outdat = \
          self.sum_reduce_all_data(list_of_local_data, t_output, comm)
        if self.verbose:
            dhwdt_abs2_locsum = np.sum(list_of_dhwdt_abs2, axis=0)
            dhwdt_abs2_totals = np.zeros_like(dhwdt_abs2_locsum)\
              if rank == root else None

        if self.sitedata:
            sij = self.sum_reduce_site_data(list_of_local_ijdata, t_output,\
              self.tpnt_sites, comm)
            if rank == root:
                sij.normalize_data(self.n_t)
                sij.dump_data()

        if self.verbose:
            temp_comm = Intracomm(comm)
            dhwdt_abs2_totals = temp_comm.reduce(dhwdt_abs2_locsum, root=root)
            if rank == root:
                dhwdt_abs2_totals = dhwdt_abs2_totals/(self.n_t * N * N)
                dhwdt_abs_totals = np.sqrt(dhwdt_abs2_totals)

        #Dump to file
        if rank == root:
            outdat.normalize_data(self.n_t, N)
            if self.file_output:
                outdat.dump_data()
            if self.verbose:
                print("t-deriv of Hamilt (abs square) with wigner avg: ")
                print("  ")
                print(tabulate({"time": t_output, \
                  "dhwdt_abs": dhwdt_abs_totals}, \
                  headers="keys", floatfmt=".6f"))
            if self.jac and self.verbose:
                print('# Cumulative number of Jacobian evaluations by root:', \
                  np.sum(info['nje']))
            print('# Done!')
            np.seterr(**old_settings)  # reset to default
            return outdat
        else:
            np.seterr(**old_settings)  # reset to default
            return None

    def evolve(self, time_info, sampling="spr"):
        if self.s_order:
            return self.dtwa_ising_longrange_2ndorder(time_info, sampling)
        else:
            return self.dtwa_ising_longrange_1storder(time_info)

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    #Initiate the parameters in object
    p = ParamData(latsize=101, beta=1.0)
    #Initiate the DTWA system with the parameters and niter
    d = Dtwa_System(p, comm, n_t=2000)
    data = d.evolve((0.0, 1.0, 1000))
