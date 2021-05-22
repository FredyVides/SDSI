#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 02:57:52 2021
SPSOLVER  Sparse linear regression solver
   Code by Fredy Vides
   For Paper, "Sparse system identification by low-rank approximation"
   by F. Vides
@author: Fredy Vides
"""

def spsolver(A,Y,L=100,tol=1e-2,delta=1e-2):
    from numpy.linalg import svd,lstsq,norm
    from numpy import zeros,dot,diag,argsort,sort,inf
    N=A.shape[1]
    M=Y.shape[1]
    X=zeros((M,N))
    u,s,v=svd(A,full_matrices=0)
    rk=sum(s>tol)
    u=u[:,:rk]
    s=s[:rk]
    s=1/s
    s=diag(s)
    v=v[:rk,:]
    A=dot(u.T,A)
    Y=dot(u.T,Y)
    X0=dot(v.T,dot(s,Y))
    w=zeros((N,))
    for k in range(M):
        K=1
        Error=1+tol
        c=X0[:,k]
        x0=c
        ac=abs(c)
        f=argsort(-ac)
        N0=max(sum(ac[f]>delta),1)
        while (K<=L) & (Error>tol):
            ff=sort(f[:N0])
            X[:,k]=w
            c, res, rnk, s = lstsq(A[:,ff],Y[:,k],rcond=None)
            X[ff,k]=c
            Error=norm(x0-X[:,k],inf)
            x0=X[:,k]
            ac=abs(x0)
            f=argsort(-ac)
            N0=max(sum(ac[f]>delta),1)
            K=K+1
    return X
