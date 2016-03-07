

# program used to generate Fig. 2, supplementary Fig. 1, and supplementari Fig. 2 of the paper
# "Free energy of ligand-receptor systems forming multimeric complexes",
# Authors: L. Di Michele, S. Bachman, L. Parolini, and B. M. Mognetti


import numpy as np

# L: length of the DNA spacers in unit of nm
L = 10.0
# R: radius of colloids in unit of L
R = 10.0
# Na: Avogadro number
Na = 6.02214179*10**23;
# rho0: standard concentration in unit of L
rho0 = Na/(10**(-3) * 10**27)*L**3;
OmegaInf=4*np.pi/3*((R+1)**3-R**3);
# Ntot: number of DNA strands of a given family on each particle (the particle of type A is functionalised by 3 Ntot strands in total)
Ntot=100.0;

# here we assume that d>|R1-R2|
def VolOvl(R1,R2,d):
    if(d>R1+R2):
        z=0
    else:
        dp=R1+R2
        dm=R1-R2
        z=np.pi/(12*d)*((dp-d)**2)*(d**2+2*d*dp-3*dm**2)
    return z

def VolA(nneigh,d):
    z=OmegaInf-nneigh*VolOvl(R+1,R,d)
    return z

def VolB(d):
    z=OmegaInf-VolOvl(R+1,R,d)
    return z

def VolAB(d):
    z=VolOvl(R+1,R+1,d)-2*VolOvl(R+1,R,d)
    return z

def FreeRep(nneigh,d):
    Vovld=VolOvl(R+1,R,d)
    VrepB=-nneigh*Ntot*np.log(1-Vovld/OmegaInf)
    VrepA=-3.*Ntot*np.log(1-nneigh*Vovld/OmegaInf)
    Vtot=VrepA+VrepB
    return Vtot

def Sigmal(nneigh,d,DGl):
    z=np.exp(-DGl)/(rho0*VolA(nneigh,d))
    return z

def Sigmas(nneigh,d,DGs):
    z=np.exp(-DGs)/((rho0*VolA(nneigh,d))**2)
    return z

def Sigmab(nneigh,d,DGb):
    z=np.exp(-DGb)*VolAB(d)/(rho0*VolA(nneigh,d)*VolB(d))
    return z

def itera(variables):
    nneigh, d, DGb, DGl, DGs = variables
    e1=1
    e2=1
    e3=1
    if(nneigh==0):
        e1=0
        e2=0
        e3=0
    if(nneigh==1):
        e2=0
        e3=0
    if(nneigh==2):
        e3=0

    db=Sigmab(nneigh,d,DGb)
    dl=Sigmal(nneigh,d,DGl)
    ds=Sigmas(nneigh,d,DGs)

    nA1o=Ntot
    nA2o=Ntot
    nA3o=Ntot
    nB1o=Ntot
    nB2o=Ntot
    nB3o=Ntot
    nA1=Ntot
    nA2=Ntot
    nA3=Ntot
    nB1=Ntot
    nB2=Ntot
    nB3=Ntot

    iit=0
    error=1
    MaxIt=10000
    MaxErr=0.0000001

    while (iit<MaxIt and error>MaxErr ):

        nA1=Ntot/(1+(nA2+nA3)*dl+nA2*nA3*ds+e1*nB1*db)
        nA2=Ntot/(1+(nA1+nA3)*dl+nA1*nA3*ds+e2*nB2*db)
        nA3=Ntot/(1+(nA1+nA2)*dl+nA1*nA2*ds+e3*nB3*db)
        nB1=Ntot/(1+e1*nA1*db)
        nB2=Ntot/(1+e2*nA2*db)
        nB3=Ntot/(1+e3*nA3*db)

        acterr=np.absolute(nA1-nA1o)/Ntot
        if(np.absolute(nA1-nA1o)/Ntot>acterr): acterr=np.absolute(nA1-nA1o)/Ntot
        if(np.absolute(nA2-nA2o)/Ntot>acterr): acterr=np.absolute(nA2-nA2o)/Ntot
        if(np.absolute(nA3-nA3o)/Ntot>acterr): acterr=np.absolute(nA3-nA3o)/Ntot
        if(np.absolute(nB1-nB1o)/Ntot>acterr): acterr=np.absolute(nB1-nB1o)/Ntot
        if(np.absolute(nB2-nB2o)/Ntot>acterr): acterr=np.absolute(nB2-nB2o)/Ntot
        if(np.absolute(nB3-nB3o)/Ntot>acterr): acterr=np.absolute(nB3-nB3o)/Ntot
        
        nA1o=nA1
        nA2o=nA2
        nA3o=nA3
        nB1o=nB1
        nB2o=nB2
        nB3o=nB3

        error=acterr
        iit=iit+1

    z = np.array([nA1, nA2, nA3, nB1, nB2, nB3])

    return z


def selfconsistent(variables):
    DGl, DGs, DGb, d, nneigh = variables
    
    z=itera([nneigh, d, DGb, DGl, DGs])
    
    nA1=z[0]
    nA2=z[1]
    nA3=z[2]
    nB1=z[3]
    nB2=z[4]
    nB3=z[5]
    
    nl1=nA2*nA3*Sigmal(nneigh,d,DGl)
    nl2=nA1*nA3*Sigmal(nneigh,d,DGl)
    nl3=nA1*nA2*Sigmal(nneigh,d,DGl)
    ns=nA1*nA2*nA3*Sigmas(nneigh,d,DGs)
    nb1=nA1*nB1*Sigmab(nneigh,d,DGb)
    nb2=nA2*nB2*Sigmab(nneigh,d,DGb)
    nb3=nA3*nB3*Sigmab(nneigh,d,DGb)
    if(nneigh==0):
        nb1=0
        nb2=0
        nb3=0
    if(nneigh==1):
        nb2=0
        nb3=0
    if(nneigh==2):
        nb3=0

    z=[nl1, nl2, nl3, ns, nb1, nb2, nb3]

    return z

def freeenergy(variables):
    nl1, nl2, nl3, ns, nb1, nb2, nb3 = variables

    if(nneigh==0):
        nb1=0
        nb2=0
        nb3=0
    if(nneigh==1):
        nb2=0
        nb3=0
    if(nneigh==2):
        nb3=0

    nA1=Ntot-nl2-nl3-ns-nb1
    nA2=Ntot-nl1-nl3-ns-nb2
    nA3=Ntot-nl1-nl2-ns-nb3

    nB1=Ntot-nb1
    nB2=Ntot-nb2
    nB3=Ntot-nb3

    free=np.log(nB1/Ntot)+np.log(nB2/Ntot)+np.log(nB3/Ntot)
    free+=np.log(nA1/Ntot)+np.log(nA2/Ntot)+np.log(nA3/Ntot)
    free*=Ntot
    free+=nl1+nl2+nl3+nb1+nb2+nb3
    free+=2*ns

    return free

########################################
# part of code used to draw supplementary Fig. 1
########################################

"""
loop1 = []
loop2 = []
loop3 = []
loop = []
spider = []
bridge = []
valDGl = []
npoint=100
for ii in range(0,npoint,1):
    DGl=-11+11.*ii/npoint
    DGs=3*DGl
    DGb=-100
    d=2*R+1
    nneigh=0
    sol=selfconsistent([DGl, DGs, DGb, d, nneigh])
    valDGl.append(DGl)
    loop1.append(sol[0])
    loop2.append(sol[1])
    loop3.append(sol[2])
    loop.append(sol[0]+sol[1]+sol[2])
    spider.append(sol[3])
    bridge.append(sol[4]+sol[5]+sol[6])

import matplotlib.pyplot as plt
plt.figure(1)
plt.xlabel('DGl')
plt.ylabel('number of costructs')
plt.plot(valDGl, loop, label='loop')
#plt.plot(valDGl, loop2, label='loop2')
#plt.plot(valDGl, loop3, label='loop3')
plt.plot(valDGl, spider, label='spider')
plt.legend()
#plt.figure(2)
#plt.plot(valDGl, bridge, label='bridge')
#plt.legend()
plt.show()
"""

######################################################
# part of code used to draw Fig. 2b and supplementary Fig. 2
######################################################
"""
bridge1 = []
bridge2 = []
bridge3 = []
spider = []
loop1 = []
loop2 = []
loop3 = []
free = []
dfree = []
valDGb =[]
nneighList = [0, 1, 2, 3];
DGl=-10.0
DGs=3*DGl
d=2*R+1

for nneigh in nneighList:

    bridge1.append([])
    bridge2.append([])
    bridge3.append([])
    spider.append([])
    loop1.append([])
    loop2.append([])
    loop3.append([])
    free.append([])

    for DGb in range(0,-30,-1):
        
        if(nneigh==0): valDGb.append(DGb)
        sol=selfconsistent([DGl, DGs, DGb, d, nneigh])
        loop1[nneigh].append(sol[0])
        loop2[nneigh].append(sol[1])
        loop3[nneigh].append(sol[2])
        spider[nneigh].append(sol[3])
        bridge1[nneigh].append(sol[4])
        bridge2[nneigh].append(sol[5])
        bridge3[nneigh].append(sol[6])
        free[nneigh].append(freeenergy(sol))



# calculate difference in free between cluster with different numbers of neighbours

for nneigh in nneighList:
    ii=0
    dfree.append([])
    for DGb in range(0,-30,-1):
        if(nneigh>0):
            dfree[nneigh].append(free[nneigh][ii]-free[nneigh-1][ii])
        else:
            dfree[nneigh].append(0)
        ii=ii+1

#ii=0
#for DGb in range(0,-30,-1):
#    print valDGb[ii], dfree[1][ii], dfree[2][ii], dfree[3][ii]
#    ii+=1

import matplotlib.pyplot as plt
plt.figure(1)
plt.xlabel('DGb')
plt.ylabel('number of bridges')
plt.plot(valDGb, bridge1[1], label='1 neighbour')
plt.plot(valDGb, bridge2[2], label='2 neighbour')
plt.plot(valDGb, bridge3[3], label='3 neighbour')
plt.legend()

plt.figure(2)
plt.xlabel('DGb')
plt.ylabel('number of spiders')
plt.plot(valDGb, spider[1], label='1 neighbour')
plt.plot(valDGb, spider[2], label='2 neighbour')
plt.plot(valDGb, spider[3], label='3 neighbour')
plt.legend(loc=4)

plt.figure(3)
plt.xlabel('DGb')
plt.ylabel('Free(i)-Free(i-1)')
plt.plot(valDGb, dfree[1], label='Free(AB)-Free(A)')
plt.plot(valDGb, dfree[2], label='Free(AB2)-Free(AB1)')
plt.plot(valDGb, dfree[3], label='Free(AB3)-Free(AB2)')
plt.legend(loc=4)
plt.show()
"""

################################################
# part of code used to draw Manuscript Fig. 2 c
################################################


free = []
dfree = []
vald = []
nneighList = [0, 1, 2, 3];
DGl=-10.0
DGs=3*DGl
DGb=-14

for nneigh in nneighList:
    free.append([])
    for id in range(1,200, 1):
        d=2*R+0.01*id
        if(nneigh==0): vald.append(d-2.*R)
        sol=selfconsistent([DGl, DGs, DGb, d, nneigh])
        free[nneigh].append(freeenergy(sol)+FreeRep(nneigh,d))

# calculate difference in free between cluster with different numbers of neighbours

for nneigh in nneighList:
    ii=0
    dfree.append([])
    for id in range(1,200, 1):
        if(nneigh>0):
            dfree[nneigh].append(free[nneigh][ii]-free[nneigh-1][ii])
        else:
            dfree[nneigh].append(0)
        ii=ii+1

for id in range(1,200, 1):
    d=vald[id-1]
    print d, dfree[1][id-1], dfree[2][id-1], dfree[3][id-1]

import matplotlib.pyplot as plt
plt.figure(1)
plt.xlabel('(d-2R)/L')
plt.ylabel(r'$\Delta F$')
plt.plot(vald, dfree[1], label='Free(AB)-Free(A)')
plt.plot(vald, dfree[2], label='Free(AB2)-Free(AB1)')
plt.plot(vald, dfree[3], label='Free(AB3)-Free(AB2)')
plt.legend(loc=4)
plt.show()

