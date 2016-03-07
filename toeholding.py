
# program used to generate Fig. 3 of the paper

# "Free energy of ligand-receptor systems forming multimeric complexes",
# Authors: L. Di Michele, S. Bachman, L. Parolini, and B. M. Mognetti



from scipy.optimize import fsolve
import numpy as np

Na = 6.02214179*10**23;
Kcal = 4186.7999409;
KJ = 1000.;
Kb = 1.3806505*10**(-23);
rho0 = Na/(10**(-3) * 10**27);

A=4.*np.pi*(200.**2)
L=10.
h=14.
conf=1./(A*L*rho0)
molt=5.

kcal = 1.74755

def Keqb1(DG1):
    z=conf*(2-h/L)*np.exp(-DG1*kcal)
    return z

def Keql1(DG1):
    z=conf*np.exp(-DG1*kcal)
    return z

def Keqb2(DG2):
    z=conf*(2-h/L)*np.exp(-DG2*kcal)
    return z

def Keql2(DG2):
    z=conf*np.exp(-DG2*kcal)
    return z

def Keqtb(DGtre):
    z=molt*(conf**2)*(2-h/L)*np.exp(-DGtre*kcal)
    return z

def Keqtl(DGtre):
    z=molt*(conf**2)*np.exp(-DGtre*kcal)
    return z

# in unit of Kcal/mol from DINAMelt without dangles
def DG1(T):
    z = -14.9251 + 0.1771 * T
    return z

def DG2(T):
    z = -14.4863 + 0.1615 * T
    return z

def DG3(T):
    z = -19.2445 + 0.24322 * T
    return z


def itera(variables):
    N, T, chi1, chi2 = variables

    nA1=N*chi1
    nA2=N*chi2
    nB=N

    nA1o=nA1
    nA2o=nA2
    nBo=nB

    iit=0
    error=1
    MaxIt=1000000000
    MaxErr=0.00000000001

    while ( iit<MaxIt and error>MaxErr ):


        nA1=N*chi1/( 1 + Keqb1(DG1(T))*nB + Keql1(DG1(T))*nB + 3 * Keqtb(DG3(T))*nA2*nB + Keqtl(DG3(T))*nA2*nB )
        nA2=N*chi2/( 1 + Keqb2(DG2(T))*nB + Keql2(DG2(T))*nB + 3 * Keqtb(DG3(T))*nA1*nB + Keqtl(DG3(T))*nA1*nB )
        nB=N/( 1 + Keqb1(DG1(T))*nA1 + Keqb2(DG2(T))*nA2 + Keql1(DG1(T))*nA1+ Keql2(DG2(T))*nA2 + 3 * Keqtb(DG3(T))*nA1*nA2 + Keqtl(DG3(T))*nA1*nA2 )


        acterr=np.absolute((nA1-nA1o)/N)
        if(np.absolute((nA2-nA2o)/N)>acterr): acterr=np.absolute((nA2-nA2o)/N)
        if(np.absolute((nB-nBo)/N)>acterr): acterr=np.absolute((nB-nBo)/N)
        
        nA1o=nA1
        nA2o=nA2
        nBo=nB
        
        error=acterr
        
        iit=iit+1

        z = np.array([nA1, nA2, nB])
    
    return z



def itera0(variables):
    N, T, chi1, chi2 = variables
    
    nA1=N*chi1
    nA2=N*chi2
    nB=N
    
    nA1o=nA1
    nA2o=nA2
    nBo=nB
    
    iit=0
    error=1
    MaxIt=1000000000
    MaxErr=0.00000000001
    
    while ( iit<MaxIt and error>MaxErr ):
        
        
        nA1=N*chi1/( 1 + Keql1(DG1(T))*nB + Keqtl(DG3(T))*nA2*nB )
        nA2=N*chi2/( 1 + Keql2(DG2(T))*nB + Keqtl(DG3(T))*nA1*nB )
        nB=N/( 1 + Keql1(DG1(T))*nA1+ Keql2(DG2(T))*nA2 + Keqtl(DG3(T))*nA1*nA2 )
        
        acterr=np.absolute((nA1-nA1o)/N)
        if(np.absolute((nA2-nA2o)/N)>acterr): acterr=np.absolute((nA2-nA2o)/N)
        if(np.absolute((nB-nBo)/N)>acterr): acterr=np.absolute((nB-nBo)/N)
        
        nA1o=nA1
        nA2o=nA2
        nBo=nB
        
        error=acterr
        
        iit=iit+1
        
        z = np.array([nA1, nA2, nB])
    
    return z



def selfconsistent(variables):

    T, N, chi1, chi2 = variables
    z=itera([N, T, chi1, chi2])
    nA1=z[0]
    nA2=z[1]
    nB=z[2]

    nb1=Keqb1(DG1(T))*nA1*nB
    nb2=Keqb2(DG2(T))*nA2*nB
    nl1=Keql1(DG1(T))*nA1*nB
    nl2=Keql2(DG2(T))*nA2*nB

    ntb=Keqtb(DG3(T))*nA1*nA2*nB
    ntl=Keqtl(DG3(T))*nA1*nA2*nB

    zz=np.array([nb1, nb2, nl1, nl2, ntb, ntl, chi1, chi2])

    return zz


def selfconsistent0(variables):
    
    T, N, chi1, chi2 = variables
    z=itera0([N, T, chi1, chi2])
    nA1=z[0]
    nA2=z[1]
    nB=z[2]
    
    nb1=0.
    nb2=0.
    nl1=Keql1(DG1(T))*nA1*nB
    nl2=Keql2(DG2(T))*nA2*nB
    
    ntb=0.
    ntl=Keqtl(DG3(T))*nA1*nA2*nB
    
    zz=np.array([nb1, nb2, nl1, nl2, ntb, ntl, chi1, chi2])
    
    return zz



def freeenergy(variables):

    b1, b2, l1, l2, tb, tl, chi1, chi2 = variables
    B=N-b1-b2-l1-l2-3*tb-tl
    A1=chi1*N-b1-l1-3*tb-tl
    A2=chi2*N-b2-l2-3*tb-tl
    z=b1+b2+l1+l2+6*tb+2*tl
    if(B>0): z=z+N*np.log(B/N)
    if(A1>0): z=z+chi1*N*np.log(A1/(N*chi1))
    if(A2>0): z=z+chi2*N*np.log(A2/(N*chi2))
    return z


N=360
Tlist=[15,20,25,30,35]
two =[]
three = []
free = []
chilist = []

for ilist in range(5):
    T=Tlist[ilist]
    two.append([])
    three.append([])
    free.append([])
    for ilist2 in range(101):
        chi1=2.*ilist2/100
        if(ilist==0): chilist.append(chi1)
        AA=selfconsistent([T,N,chi1,2-chi1])
        AA0=selfconsistent0([T,N,chi1,2-chi1])
        two[ilist].append((AA[0]+AA[1]+AA[2]+AA[3])/N)
        three[ilist].append((3*AA[4]+AA[5])/N)
        free[ilist].append(freeenergy(AA)/N-freeenergy(AA0)/N)


import matplotlib.pyplot as plt
plt.figure(1)
plt.xlabel('chi')
plt.ylabel('two')
plt.plot(chilist, two[0], label='T=15')
plt.plot(chilist, two[1], label='T=20')
plt.plot(chilist, two[2], label='T=25')
plt.plot(chilist, two[3], label='T=30')
plt.plot(chilist, two[4], label='T=35')
plt.legend()
plt.figure(2)
plt.xlabel('chi')
plt.ylabel('three')
plt.plot(chilist, three[0], label='T=15')
plt.plot(chilist, three[1], label='T=20')
plt.plot(chilist, three[2], label='T=25')
plt.plot(chilist, three[3], label='T=30')
plt.plot(chilist, three[4], label='T=35')
plt.legend()
plt.figure(3)
plt.xlabel('chi')
plt.ylabel('free')
plt.plot(chilist, free[0], label='T=15')
plt.plot(chilist, free[1], label='T=20')
plt.plot(chilist, free[2], label='T=25')
plt.plot(chilist, free[3], label='T=30')
plt.plot(chilist, free[4], label='T=35')
plt.legend()
plt.show()



