#!/usr/bin/env python
"""
killMS, a package for calibration in radio interferometry.
Copyright (C) 2013-2017  Cyril Tasse, l'Observatoire de Paris,
SKA South Africa, Rhodes University

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import optparse
import sys
from killMS.Other import MyPickle
from killMS.Other import logo
from killMS.Other import ModColor
from DDFacet.Other import logger
log=logger.getLogger("killMS")
#logger.itsLog.logger.setLevel(logger.logging.CRITICAL)

sys.path=[name for name in sys.path if not(("pyrap" in name)&("/usr/local/lib/" in name))]

# test

#import numpy
#print numpy.__file__
#import pyrap
#print pyrap.__file__
#stop


if "nocol" in sys.argv:
    print("nocol")
    ModColor.silent=1
    
if "nox" in sys.argv:
    import matplotlib
    matplotlib.use('agg')
    print(ModColor.Str(" == !NOX! =="))

import time
import os
import numpy as np
import pickle
from itertools import product as ItP

NameSave="last_plotSols.obj"
def read_options():

    desc="""killMS Questions and suggestions: cyril.tasse@obspm.fr"""
    
    opt = optparse.OptionParser(usage='Usage: %prog --ms=somename.MS <options>',version='%prog version 1.0',description=desc)
    group = optparse.OptionGroup(opt, "* Data-related options", "Won't work if not specified.")
    group.add_option('--SolsFile',help='Input Solutions list [no default]',default='')
    group.add_option('--DoResid',type="int",help='No [no default]',default=-1)
    group.add_option('--PlotMode',type='str',help=' [no default]',default="AP")
    group.add_option('--DirList',type="str",help=' [no default]',default="")
    group.add_option('--NoTicks',action="store_true",help="disable Y axis scale")
    group.add_option('--ChanList',help=' [no default]',default="")
    opt.add_option_group(group)
    
    options, arguments = opt.parse_args()
    f = open(NameSave,"wb")
    pickle.dump(options,f)
    
import pylab

import numpy as np

def GiveNXNYPanels(Ns,ratio=800/500):
    nx=int(round(np.sqrt(Ns/ratio)))
    ny=int(nx*ratio)
    if nx*ny<Ns: ny+=1
    return nx,ny

from killMS.Array import ModLinAlg

def NormMatrices(G):
    #print "no norm"
    #return G
    nt,nch,na,_,_=G.shape

    for iChan,it in ItP(range(nch),range(nt)):
        Gt=G[it,iChan,:,:]
        u,s,v=np.linalg.svd(Gt[0])
        # #J0/=np.linalg.det(J0)
        # J0=Gt[0]
        # JJ=np.dot(J0.T.conj(),J0)
        # sqJJ=ModLinAlg.sqrtSVD(JJ)
        # sqJJinv=ModLinAlg.invSVD(JJ)
        # U=np.dot(J0,sqJJinv)
        U=np.dot(u,v)
        for iAnt in range(0,na):
            Gt[iAnt,:,:]=np.dot(U.T.conj(),Gt[iAnt,:,:])
            #Gt[iAnt,:,:]=np.dot(np.dot(u,Gt[iAnt,:,:]),v.T.conj())
            #Gt[iAnt,:,:]=np.dot(Gt[iAnt,:,:],J0)
    return G


def main(options=None):
    

    if options==None:
        f = open(NameSave,'rb')
        options = pickle.load(f)



    
    FilesList=options.SolsFile.split(",")


    


    LSols=[]
    nSol=len(FilesList)
    t0=None
    for FileName in FilesList:
        
        if "npz" in FileName:
            SolsDico=np.load(FileName)
            Sols=SolsDico["Sols"]
            StationNames=SolsDico["StationNames"]
            ClusterCat=SolsDico["ClusterCat"]
            Sols=Sols.view(np.recarray)
            nt,nch,na,nd,_,_=Sols.G.shape

        elif "h5" in FileName:
            import tables
            H5=tables.openFile(FileName)
            npol, nch, nd, na, nchan, nt=H5.root.sol000.amplitude000.val.shape
            GH5=H5.root.sol000.amplitude000.val[:]*np.exp(1j*H5.root.sol000.phase000.val[:])
            Times=H5.root.sol000.amplitude000.time[:]
            StationNames=H5.root.sol000.antenna[:]["name"]
            H5.close()
            

            Sols=np.zeros((nt,),dtype=[("t0",np.float64),
                                       ("t1",np.float64),
                                       ("G",np.complex64,(nch,na,nd,2,2))])
            Sols=Sols.view(np.recarray)
            dt=np.median(Times[1::]-Times[0:-1])
            Sols.t0=Times-dt/2.
            Sols.t1=Times+dt/2.
            for iTime in range(nt):
                for iDir0,iDir1 in zip(range(3),range(3)):#[0,2,1]):#range(nd):
                    for iAnt in range(na):
                        for ipol in range(4):
                            Sols.G[iTime,iAnt,iDir0].flat[ipol]=GH5[ipol,iDir1,iAnt,0,iTime]

            
        ind=np.where(Sols.t1!=0)[0]
        Sols=Sols[ind]
        tm=(Sols.t1+Sols.t0)/2.
        if t0==None:
            t0=tm[0]
        tm-=t0
        Sols.t0=tm
        nx,ny=GiveNXNYPanels(na)
        LSols.append(Sols)

        # LSols=[LSols[0]]
        # nSol=1

    # diag terms
    Lls=["-",":",":"]
    Lcol0=["black","black","blue"]
    Lcol1=["gray","gray","red"]
    Lalpha0=[1,1,1]
    Lalpha1=[0.5,0.5,0.5]

    # Lls=["-","--",":"]
    # Lcol0=["red","blue","blue"]
    # Lcol1=["gray","red","red"]
    # # Lalpha0=[1,0.5,1]
    # # Lalpha1=[0.5,0.5,0.5]


    # off-diag terms
    Lls_off=Lls#["-","--",":"]
    Lcol0_off=Lcol0#["black","black","blue"]
    Lcol1_off=Lcol1#["gray","gray","red"]
    
    if options.DoResid!=-1:
        Sresid=LSols[1].copy()
        LSols.append(Sresid)

    if options.DirList!="":
        DirList=options.DirList.split(',')
        DirList=[int(i) for i in DirList]
    else:
        DirList=range(nd)

    if options.ChanList!="":
        ChanList=options.ChanList.split(',')
        ChanList=[int(i) for i in ChanList]
    else:
        ChanList=range(nch)

    #DirList=[np.where(ClusterCat["SumI"]==np.max(ClusterCat["SumI"]))[0][0]]
    #print DirList

    for iDir in DirList:
        pylab.figure(0,figsize=(13,8))
        iAnt=0
        for iSol in range(nSol):
            Sols=LSols[iSol]
            G=Sols.G[:,:,:,iDir,:,:]
            Sols.G[:,:,:,iDir,:,:]=NormMatrices(G)
            
        ampMax=1.5*np.max(np.median(np.abs(LSols[0].G),axis=0))
        if options.PlotMode=="AP":
            op0=np.abs
            op1=np.angle
            ylim0=0,ampMax
            ylim1=-np.pi,np.pi
            PlotDiag=[True,False]
        elif options.PlotMode=="ReIm":
            op0=np.real
            op1=np.imag
            ylim0=-ampMax,ampMax
            ylim1=-ampMax,ampMax
            PlotDiag=[True,True]
        elif options.PlotMode=="A":
            op0=np.abs
            op1=None
            ylim0=0,ampMax
            PlotDiag=[True]
        elif options.PlotMode=="P":
            op0=np.angle
            op1=None
            ylim0=-np.pi,np.pi
            PlotDiag=[False]

        L_ylim0=(0,1.5*np.max(np.median(np.abs(LSols[0].G[:,:,:,iDir,:,:]),axis=0)))

        if options.DoResid!=-1:
            LSols[-1].G[:,:,iDir,:,:]=LSols[1].G[:,:,iDir,:,:]-LSols[0].G[:,:,iDir,:,:]
            nSol+=1
            
        marker="."
        marker=""

        pylab.clf()
        for i in range(nx):
            for j in range(ny):
                if iAnt>=na:continue
                if iAnt>=1:
                    ax=pylab.subplot(nx,ny,iAnt+1,sharex=axRef,sharey=axRef)
                else:
                    axRef=pylab.subplot(nx,ny,iAnt+1)
                    ax=axRef

                ax.set_xticks([])
                if options.NoTicks:
                    ax.set_yticks([])

                ax.tick_params(axis='both', which='major', labelsize='xx-small')
                if options.PlotMode == "P":
                    ax.set_ylim(ylim0)
                else:
                    ax.set_ylim(L_ylim0)

                if op1 is not None:
                    ax2 = ax.twinx()
                    ax2.set_ylim(ylim1)
                    ax2.set_xticks([])
                    if options.NoTicks or iAnt % ny != ny - 1:
                        ax2.set_yticks([])
                    ax2.tick_params(axis='both', which='major', labelsize='xx-small')

                pylab.title(StationNames[iAnt], fontsize=9)
                for iSol in range(nSol):
                    Sols=LSols[iSol]
                    ThisChanList=range(Sols.G.shape[1])
                    for iChan in ThisChanList:                    
                        G=Sols.G[:,iChan,:,iDir,:,:]
                        J=G[:,iAnt,:,:]
                        if op1!=None:
                            # ax.plot(tm,op1(J[:,0,1]),color="blue")
                            # ax.plot(tm,op1(J[:,1,0]),color="blue")
                            ax2.plot(Sols.t0,op1(J[:,1,1]),color=Lcol1[iSol],alpha=Lalpha1[iSol],ls=Lls[iSol],marker=marker,zorder=0)
                            ax2.plot(Sols.t0,op1(J[:,0,0]),color=Lcol1[iSol],alpha=Lalpha1[iSol],ls=Lls[iSol],marker=marker,zorder=0)
                            if PlotDiag[1]:
                                ax2.plot(Sols.t0,op1(J[:,0,1]),color=Lcol1_off[iSol],alpha=Lalpha1[iSol],ls=Lls_off[iSol],marker=marker,zorder=0)
                                ax2.plot(Sols.t0,op1(J[:,1,0]),color=Lcol1_off[iSol],alpha=Lalpha1[iSol],ls=Lls_off[iSol],marker=marker,zorder=0)
                            #print StationNames[iAnt]
                        
                        ax.plot(Sols.t0,op0(J[:,0,0]),color=Lcol0[iSol],alpha=Lalpha0[iSol],ls=Lls[iSol],marker=marker,zorder=10)
                        ax.plot(Sols.t0,op0(J[:,1,1]),color=Lcol0_off[iSol],alpha=Lalpha0[iSol],ls=Lls_off[iSol],marker=marker,zorder=10)
                        if PlotDiag[0]:
                            ax.plot(Sols.t0,op0(J[:,1,0]),color=Lcol0[iSol],alpha=Lalpha0[iSol],ls=Lls[iSol],marker=marker,zorder=10)
                            ax.plot(Sols.t0,op0(J[:,0,1]),color=Lcol0_off[iSol],alpha=Lalpha0[iSol],ls=Lls_off[iSol],marker=marker,zorder=10)
                        if options.PlotMode=="P":
                            ax.set_ylim(ylim0)
                        else:
                            ax.set_ylim(L_ylim0)
    
                        ax.set_xticks([])
                        if options.NoTicks:
                            ax.set_yticks([])


                iAnt+=1
        pylab.suptitle('Direction %i [%s]'%(iDir,str(ClusterCat["SumI"][iDir])))#L_ylim0)))
        pylab.tight_layout(pad=3., w_pad=0.5, h_pad=2.0)
        pylab.draw()
        pylab.show()
        #pylab.pause(0.1)
        #time.sleep(1)


if __name__=="__main__":
    read_options()
    f = open(NameSave,'rb')
    options = pickle.load(f)

    main(options=options)
