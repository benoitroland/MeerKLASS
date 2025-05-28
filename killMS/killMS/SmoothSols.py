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
import pickle
import numpy as np
import numpy as np
#import pylab
import os
from DDFacet.Other import logger
from DDFacet.Other import ModColor
log=logger.getLogger("ClassInterpol")
from DDFacet.Other import AsyncProcessPool
from DDFacet.Other import Multiprocessing
from DDFacet.Array import shared_dict
from killMS.Array import NpShared
IdSharedMem=str(int(os.getpid()))+"."
from DDFacet.Other import AsyncProcessPool
from killMS.Other import ClassFitTEC
#from killMS.Other import ClassFit_1D
from killMS.Other import ClassFitAmp
#from killMS.Other import ClassFitAmp
from killMS.Other import ClassClip
import scipy.ndimage.filters
from pyrap.tables import table
# # ##############################
# # Catch numpy warning
# np.seterr(all='raise')
# import warnings
# #with warnings.catch_warnings():
# #    warnings.filterwarnings('error')
# warnings.catch_warnings()
# warnings.filterwarnings('error')
# # ##############################
from killMS.Other.ClassTimeIt import ClassTimeIt
#from killMS.Other.least_squares import least_squares
from scipy.optimize import least_squares
import copy
from pyrap.tables import table
SaveName="last_InterPol.obj"
from itertools import product as ItP

def NormMatrices(G):
    log.print("Normalising to a reference antenna...")
    nt, nch, na, nDir, _, _=G.shape
    for iDir in range(nDir):
        Gg=G[:,:,:,iDir,:,:]
        G[:,:,:,iDir,:,:]=NormMatricesDir(Gg)

def NormMatricesDir(G):
    nt,nch,na,_,_=G.shape
    for iChan,it in ItP(range(nch),range(nt)):
        Gt=G[it,iChan,:,:]
        u,s,v=np.linalg.svd(Gt[0])
        U=np.dot(u,v)
        for iAnt in range(0,na):
            Gt[iAnt,:,:]=np.dot(U.T.conj(),Gt[iAnt,:,:])
    return G

def read_options():
    desc="""Questions and suggestions: cyril.tasse@obspm.fr"""
    global options
    opt = optparse.OptionParser(usage='Usage: %prog --ms=somename.MS <options>',version='%prog version 1.0',description=desc)

    group = optparse.OptionGroup(opt, "* Data-related options", "Won't work if not specified.")
    group.add_option('--SolsFileIn',help='SolfileIn [no default]',default=None)
    group.add_option('--SolsFileOut',help='SolfileOut [no default]',default=None)
    group.add_option('--InterpMode',help='Interpolation mode TEC and/or Amp [default is %default]',type="str",default="TEC,Amp")
    group.add_option('--CrossMode',help='Use cross gains maode for TEC [default is %default]',type=int,default=1)
    group.add_option('--RemoveAmpBias',help='Remove amplitude bias (along time) before smoothing [default is %default]',type=int,default=0)
    group.add_option('--RemoveMedianAmp',help='Remove median amplitude (along freq) after fitting [default is %default]',type=int,default=1)
    
    group.add_option('--Amp-SmoothType',help='Interpolation Type for the amplitude [default is %default]',type="str",default="Gauss")
    group.add_option('--Amp-PolyOrder',help='Order of the polynomial to do the amplitude',type="int",default=2)
    group.add_option('--Amp-GaussKernel',help='',type="str",default="1,3")
    group.add_option('--NCPU',help='Number of CPU to use',type="int",default=0)
    opt.add_option_group(group)


    options, arguments = opt.parse_args()
    exec("options.Amp_GaussKernel=%s"%options.Amp_GaussKernel)
    f = open(SaveName,"wb")
    pickle.dump(options,f)

def TECToPhase(TEC,freq):
    K=8.4479745e9
    phase=K*TEC*(1./freq)
    return phase

def TECToZ(TEC,ConstPhase,freq):
    return np.exp(1j*(TECToPhase(TEC,freq)+ConstPhase))

class ClassInterpol():
    def __init__(self,
                 InSolsName=None,
                 OutSolsName=None,
                 InterpMode="TEC",
                 PolMode="Scalar",
                 Amp_PolyOrder=1,NCPU=0,
                 Amp_GaussKernel=(0,5),
                 Amp_SmoothType="Poly",
                 CrossMode=1,
                 RemoveAmpBias=0,
                 DoPBar=True,
                 RemoveMedianAmp=True):

        
        if type(InterpMode)==str:
            InterpMode=InterpMode.split(",")#[InterpMode]
        self.DoPBar=DoPBar
        self.InSolsName=InSolsName
        self.RemoveAmpBias=RemoveAmpBias
        if self.InSolsName is not None:
            self.setInSols(InSolsName)
            
        self.OutSolsName=OutSolsName
        self.RemoveMedianAmp=RemoveMedianAmp
        self.workersHaveStarted=False
        

        self.CrossMode=CrossMode
        self.incrCross=11


        self.InterpMode=InterpMode
        self.Amp_PolyOrder=Amp_PolyOrder


        
        self.PolMode=PolMode
        self.Amp_GaussKernel=Amp_GaussKernel
        if len(self.Amp_GaussKernel)!=2:
            raise ValueError("GaussKernel should be of size 2")
        self.Amp_SmoothType=Amp_SmoothType

        self.ModeFitTEC=["TEC","CPhase"]
        if "TEC" in self.InterpMode:
            if "TEC" in self.ModeFitTEC:
                log.print( "  Smooth phases using a TEC model")
            if "CPhase" in self.ModeFitTEC:
                log.print( "  Smooth phases using a Constant Phase model")
            if self.CrossMode: 
                log.print(ModColor.Str("Using CrossMode"))
            

                
        if "Amp" in self.InterpMode:
            if Amp_SmoothType=="Poly":
                log.print( "  Smooth amplitudes using polynomial model of order %i"%self.Amp_PolyOrder)
            if Amp_SmoothType=="Gauss":
                log.print( "  Smooth amplitudes using Gaussian kernel of %s (Time/Freq) bins"%str(Amp_GaussKernel))


        # AsyncProcessPool.init(ncpu=NCPU,affinity=0)
        # AsyncProcessPool.APP.registerJobHandlers(self)

    def setInSols(self,InSolsName):
        self.InSolsName=InSolsName
        if isinstance(self.InSolsName,str):
            log.print("Loading %s"%self.InSolsName)
            self.DicoFile=dict(np.load(self.InSolsName,allow_pickle=True))
        else:
            self.DicoFile=self.InSolsName
        
        self.Sols=self.DicoFile["Sols"].view(np.recarray)
        if "MaskedSols" in self.DicoFile.keys():
            MaskFreq=np.logical_not(np.all(np.all(np.all(self.DicoFile["MaskedSols"][...,0,0],axis=0),axis=1),axis=1))
            nt,_,na,nd,_,_=self.Sols.G.shape

            self.DicoFile["FreqDomains"]=self.DicoFile["FreqDomains"][MaskFreq]
            NFreqsOut=np.count_nonzero(MaskFreq)
            log.print("There are %i non-zero freq channels"%NFreqsOut)
            SolsOut=np.zeros((nt,),dtype=[("t0",np.float64),("t1",np.float64),
                                          ("G",np.complex64,(NFreqsOut,na,nd,2,2)),
                                          ("Stats",np.float32,(NFreqsOut,na,4))])
            SolsOut=SolsOut.view(np.recarray)
            SolsOut.G=self.Sols.G[:,MaskFreq,...]
            SolsOut.t0=self.Sols.t0
            SolsOut.t1=self.Sols.t1
            self.Sols=self.DicoFile["Sols"]=SolsOut
            del(self.DicoFile["MaskedSols"])

        MeanDataFreqPerSolBand=self.DicoFile.get("MeanDataFreqPerSolBand",None)
        if MeanDataFreqPerSolBand is not None:
            log.print("Getting frequencies for the fit at the mean freq location of the data")
            self.CentralFreqs=MeanDataFreqPerSolBand
        else:
            log.print("Getting frequencies for the fit at the center of the solution bins")
            self.CentralFreqs=np.mean(self.DicoFile["FreqDomains"],axis=1)
            
        g00=self.Sols.G[...,0,0]
        g11=self.Sols.G[...,1,1]
        g01=self.Sols.G[...,0,1]
        g10=self.Sols.G[...,1,0]
        if np.all(g00==g11) and np.alltrue(g01==0) and  np.alltrue(g10==0):
            log.print("Jones matrices are scalar...")
            self.JonesMode="Scalar"
            self.NPol=1
            self.polMap=[(0,0)]
        elif not np.all(g00==g11) and np.alltrue(g01==0) and np.alltrue(g10==0):
            log.print("Jones matrices are diagonal...")
            self.JonesMode="Diag"
            self.NPol=2
            self.polMap=[(0,0),(1,1)]
        else:
            log.print("Jones matrices are full pol...")
            self.JonesMode="Full"
            NormMatrices(self.Sols.G)
            self.NPol=4
            stop

        NTEC=101
        NConstPhase=51
        TECGridAmp=0.1
        TECGrid,CPhase=np.mgrid[-TECGridAmp:TECGridAmp:NTEC*1j,-np.pi:np.pi:NConstPhase*1j]
        Z=TECToZ(TECGrid.reshape((-1,1)),CPhase.reshape((-1,1)),self.CentralFreqs.reshape((1,-1)))
        self.Z=Z
        self.TECGrid=TECGrid
        self.CPhase=CPhase
        if self.RemoveAmpBias:
            self.G0=np.zeros((1,nch,na,nd,2,2),np.float32)
            self.CalcFreqAmpSystematics()
            self.G0[self.G0==0]=1
            self.Sols.G/=self.G0

        self.DicoStuff={}
        self.DicoStuff["GIn"]=self.Sols.G.copy()
        self.DicoStuff["GOut"]=self.Sols.G.copy()
        self.DicoStuff["Z"]=self.Z
        self.DicoStuff["TECGrid"]=self.TECGrid
        self.DicoStuff["CPhase"]=self.CPhase
        self.DicoStuff["CentralFreqs"]=self.CentralFreqs
        #self.DicoStuff["Sols"]=self.Sols
        
        if self.RemoveAmpBias:
            self.GOut*=self.G0
            
        self.DicoStuff=shared_dict.dict_to_shm("DicoStuff", self.DicoStuff)
        self.GOut=self.DicoStuff["GOut"]
        self.GIn=self.DicoStuff["GIn"]
        
        
        #self.Sols=self.Sols[0:10].copy()

    def _loadDicoStuff(self):
        self.DicoStuff=shared_dict.attach("DicoStuff")
        self.GOut=self.DicoStuff["GOut"]
        self.GIn=self.DicoStuff["GIn"]
        self.Z=self.DicoStuff["Z"]
        self.TECGrid=self.DicoStuff["TECGrid"]
        self.CPhase=self.DicoStuff["CPhase"]
        self.CentralFreqs=self.DicoStuff["CentralFreqs"]
        # self.Sols=self.DicoStuff["Sols"]
        # self.Sols=self.Sols.view(np.recarray)
        
        
    def updateSols(self,DicoSols):
        self.GIn[:]=DicoSols["Sols"]["G"][:]
        self.GOut[:]=DicoSols["Sols"]["G"][:]

        
    def TECInterPol(self):
        #Sols0=self.Sols
        nt,nch,na,nd,_,_=self.GIn.shape
        for iAnt in range(na):
            for iDir in range(nd):
                for it in range(nt):
                    self.FitThisTEC(it,iAnt,iDir)

    def CalcFreqAmpSystematics(self):
        log.print( "  Calculating amplitude systematics...")
        #Sols0=self.Sols
        nt,nch,na,nd,_,_=self.GIn.shape
        
        for iAnt in range(na):
            for iDir in range(nd):
                G=self.GIn[:,:,iAnt,iDir,:,:]
                G0=np.mean(np.abs(G),axis=0)
                self.G0[0,:,iAnt,iDir,:,:]=G0.reshape((nch,2,2))

    def startWorkers(self):
        try:
            if not self.workersHaveStarted:
                log.print("Starting workers...")
                AsyncProcessPool.APP.startWorkers()
                AsyncProcessPool.APP.awaitWorkerStart()
                self.workersHaveStarted=True
                log.print("  ... ok")
        except:
            pass


    def stopWorkers(self):
        #AsyncProcessPool.APP.terminate()
        AsyncProcessPool.APP.shutdown()
        Multiprocessing.cleanupShm()
            
    def InterpolParallel(self):
        nt,nch,na,nd,_,_=self.GIn.shape
        log.print(" #Times:      %i"%nt)
        log.print(" #Channels:   %i"%nch)
        log.print(" #Antennas:   %i"%na)
        log.print(" #Directions: %i"%nd)

        # AsyncProcessPool.APP.terminate()
        # AsyncProcessPool.APP.shutdown()
        # Multiprocessing.cleanupShm()


        
        self.startWorkers()
        
        iJob=0
        #        for iAnt in [49]:#range(na):
        #            for iDir in [0]:#range(nd):

        if "TEC" in self.InterpMode:
            #AsyncProcessPool.APP.runJob("FitThisTEC_%d"%iJob, self.FitThisTEC, args=(208,)); iJob+=1
            self.TECArray=NpShared.ToShared("%sTECArray"%IdSharedMem,np.zeros((self.NPol,nt,nd,na),np.float32))
            self.CPhaseArray=NpShared.ToShared("%sCPhaseArray"%IdSharedMem,np.zeros((self.NPol,nt,nd,na),np.float32))
            for it in range(nt):
                for iDir in range(nd):
                    AsyncProcessPool.APP.runJob("FitThisTEC_%d"%iJob, self.FitThisTECTimeDir_out, args=(it,iDir))#,serial=True)
                    iJob+=1
            progress=False
            if self.DoPBar: progress="Fit TEC"
            workers_res=AsyncProcessPool.APP.awaitJobResults("FitThisTEC*", progress=progress)

#         if "TEC" in self.InterpMode:
#             #AsyncProcessPool.APP.runJob("FitThisTEC_%d"%iJob, self.FitThisTEC, args=(208,)); iJob+=1
#             self.TECArray=NpShared.ToShared("%sTECArray"%IdSharedMem,np.zeros((nt,nd,na),np.float32))
#             self.CPhaseArray=NpShared.ToShared("%sCPhaseArray"%IdSharedMem,np.zeros((nt,nd,na),np.float32))
#             for it in range(nt):
# #            for iDir in range(nd):
#                 AsyncProcessPool.APP.runJob("FitThisTEC_%d"%iJob, self.FitThisTECTime, args=(it,))#,serial=True)
#                 iJob+=1
#             progress=False
#             if self.DoPBar: progress="Fit TEC"
#             workers_res=AsyncProcessPool.APP.awaitJobResults("FitThisTEC*", progress=progress)


        if "Amp" in self.InterpMode:
            for iAnt in range(na):
                for iDir in range(nd):
                    AsyncProcessPool.APP.runJob("FitThisAmp_%d"%iJob, self.FitThisAmp, args=(iAnt,iDir))#,serial=True)
                    iJob+=1
            progress=False
            if self.DoPBar: progress="Fit Amp"
            workers_res=AsyncProcessPool.APP.awaitJobResults("FitThisAmp*", progress=progress)

        if "PolyAmp" in self.InterpMode:
            for iDir in range(nd):
                AsyncProcessPool.APP.runJob("FitThisPolyAmp_%d"%iJob, self.FitThisPolyAmp, args=(iDir,))#,serial=True)
                iJob+=1
            progress=False
            if self.DoPBar: progress="Fit Amp"
            workers_res=AsyncProcessPool.APP.awaitJobResults("FitThisPolyAmp*", progress=progress)


        if "UnitAmp" in self.InterpMode:
            log.print("Normalising the amplitudes to 1")
            GOut=self.GOut#NpShared.GiveArray("%sGOut"%IdSharedMem)
            nt,nch,na,nd,_,_=self.GIn.shape
            for iDir in range(nd):
                for it in range(nt):
                    GOut[it,:,:,iDir,0,0]/=np.abs(GOut[it,:,:,iDir,0,0])
                    GOut[it,:,:,iDir,1,1]/=np.abs(GOut[it,:,:,iDir,1,1])

            
        if "Clip" in self.InterpMode:
            for iDir in range(nd):
                #self.ClipThisDir(iDir)
                AsyncProcessPool.APP.runJob("ClipThisDir_%d"%iJob, self.ClipThisDir, args=(iDir,))#,serial=True)
                iJob+=1
            progress=False
            if self.DoPBar: progress="Clip Amp"
            workers_res=AsyncProcessPool.APP.awaitJobResults("ClipThisDir*", progress=progress)


        # ###########################
        # import pylab
        # op0=np.abs
        # op1=np.angle
        # #for iDir in range(nd):
        # for iAnt in range(40,na):
        #     pylab.clf()
        #     A=op0(self.GIn[:,:,iAnt,iDir,0,0])
        #     v0,v1=0,A.max()
        #     pylab.subplot(2,3,1)
        #     pylab.imshow(op0(self.GIn[:,:,iAnt,iDir,0,0]),interpolation="nearest",aspect="auto",vmin=v0,vmax=v1)
        #     pylab.title("Raw Solution (Amp)")
        #     pylab.xlabel("Freq bin")
        #     pylab.ylabel("Time bin")

        #     pylab.subplot(2,3,2)
        #     pylab.imshow(op0(self.GOut[:,:,iAnt,iDir,0,0]),interpolation="nearest",aspect="auto",vmin=v0,vmax=v1)
        #     pylab.title("Smoothed Solution (Amp)")
        #     pylab.xlabel("Freq bin")
        #     pylab.ylabel("Time bin")

        #     pylab.subplot(2,3,3)
        #     pylab.imshow(op0(self.GIn[:,:,iAnt,iDir,0,0])-op0(self.GOut[:,:,iAnt,iDir,0,0]),interpolation="nearest",
        #                  aspect="auto",vmin=v0,vmax=v1)
        #     pylab.xlabel("Freq bin")
        #     pylab.ylabel("Time bin")
        #     pylab.title("Residual (Amp)")
        #     #pylab.colorbar()
        #     A=op1(self.GIn[:,:,iAnt,iDir,0,0])
        #     v0,v1=A.min(),A.max()
        #     pylab.subplot(2,3,4)
        #     pylab.imshow(op1(self.GIn[:,:,iAnt,iDir,0,0]),interpolation="nearest",aspect="auto",vmin=v0,vmax=v1)
        #     pylab.title("Raw Solution (Phase)")
        #     pylab.xlabel("Freq bin")
        #     pylab.ylabel("Time bin")

        #     pylab.subplot(2,3,5)
        #     pylab.imshow(op1(self.GOut[:,:,iAnt,iDir,0,0]),interpolation="nearest",aspect="auto",vmin=v0,vmax=v1)
        #     pylab.title("Smoothed Solution (Phase)")
        #     pylab.xlabel("Freq bin")
        #     pylab.ylabel("Time bin")

        #     pylab.subplot(2,3,6)
        #     pylab.imshow(op1(self.GIn[:,:,iAnt,iDir,0,0])-op1(self.GOut[:,:,iAnt,iDir,0,0]),
        #                  interpolation="nearest",aspect="auto",vmin=v0,vmax=v1)
        #     pylab.title("Residual (Phase)")
        #     pylab.xlabel("Freq bin")
        #     pylab.ylabel("Time bin")

        #     pylab.suptitle("(iAnt, iDir) = (%i, %i)"%(iAnt,iDir))
        #     pylab.tight_layout()
        #     pylab.draw()
        #     pylab.show()#False)
        #     pylab.pause(0.1)
        #     #stop

    
#     def FitThisTECTime(self,it):
#         self._loadDicoStuff()
#         nt,nch,na,nd,_,_=self.GIn.shape
#         TECArray=NpShared.GiveArray("%sTECArray"%IdSharedMem)
#         CPhaseArray=NpShared.GiveArray("%sCPhaseArray"%IdSharedMem)
#         for iDir in range(nd):
# #        for it in range(nt):
#             Est=None
#             if it>0:
#                 E_TEC=TECArray[it-1,iDir,:]
#                 E_CPhase=CPhaseArray[it-1,iDir,:]
#                 Est=(E_TEC,E_CPhase)
#             gz,TEC,CPhase=self.FitThisTECTimeDir(it,iDir,Est=Est)
            
#             GOut=self.GOut#NpShared.GiveArray("%sGOut"%IdSharedMem)
#             GOut[it,:,:,iDir,0,0]=gz
#             GOut[it,:,:,iDir,1,1]=gz

#             TECArray[it,iDir,:]=TEC
#             CPhaseArray[it,iDir,:]=CPhase
#             #print(it,iDir,TEC,CPhase)
#             #iAnt=-1
#             #print(TECArray[0,iDir,iAnt],GOut[0,:,iAnt,iDir,0,0])

    def FitThisTECTimeDir_out(self,it,iDir):
        self._loadDicoStuff()
        nt,nch,na,nd,_,_=self.GIn.shape
        TECArray=NpShared.GiveArray("%sTECArray"%IdSharedMem)
        CPhaseArray=NpShared.GiveArray("%sCPhaseArray"%IdSharedMem)
        
        Est=None
        polMap=self.polMap
        for iPol in range(self.NPol):
            if it>0:
                E_TEC=TECArray[iPol,it-1,iDir,:]
                E_CPhase=CPhaseArray[iPol,it-1,iDir,:]
                Est=(E_TEC,E_CPhase)
            xPol,yPol=polMap[iPol]
            gz,TEC,CPhase=self.FitThisTECTimeDirExternal(it,iDir,Est=Est,
                                                             xPol=xPol,yPol=yPol)
            GOut=self.GOut
            GOut[it,:,:,iDir,xPol,yPol]=gz
            if self.JonesMode=="Scalar":
                GOut[it,:,:,iDir,1,1]=gz
                    
            TECArray[iPol,it,iDir,:]=TEC
            CPhaseArray[iPol,it,iDir,:]=CPhase
            
            
    def FitThisTECTimeDirExternal(self,it,iDir,Est=None,xPol=0,yPol=0):
        
        GOut=self.GOut#NpShared.GiveArray("%sGOut"%IdSharedMem)
        nt,nch,na,nd,_,_=self.GIn.shape
        T=ClassTimeIt("CrossFit")
        T.disable()

        Mode=self.ModeFitTEC

        TEC0CPhase0=np.zeros((len(Mode),na),np.float32)
        for iAnt in range(na):
            _,t0,c0=self.EstimateThisTECTime(it,iAnt,iDir,xPol=xPol,yPol=yPol)
            TEC0CPhase0[0,iAnt]=t0#+np.random.randn(1)[0]*0.001
            if "CPhase" in Mode:
                TEC0CPhase0[1,iAnt]=c0

        def X2G(X):
            if "CPhase" in Mode:
                TEC,CPhase=X.reshape((len(Mode),na))
            else:
                TEC,=X.reshape((len(Mode),na))
                CPhase=np.zeros((1,na),np.float32)
            TEC-=TEC[0]
            CPhase-=CPhase[0]
            return np.abs(GOut[it,:,:,iDir,xPol,yPol]).T*TECToZ(TEC.reshape((-1,1)),CPhase.reshape((-1,1)),self.CentralFreqs.reshape((1,-1)))
        
        #GThis0=np.abs(self.GIn[it,:,:,iDir,0,0]).T.copy()#X2G(TEC0CPhase0)
        T.timeit("init")
        # ######################################
        # Changing method
        #print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        #print it,iDir
        
        if self.CrossMode:
            TECMachine=ClassFitTEC.ClassFitTEC_Cross(self.GIn[it,:,:,iDir,xPol,yPol],
                                                     self.CentralFreqs,
                                                     Tol=5.e-2,
                                                     Mode=Mode)
            TECMachine.setX0(TEC0CPhase0.ravel())
            X=TECMachine.doFit()
        else:
            LX=[]
            for iAnt in range(na):
                TECMachine=ClassFitTEC.ClassFitTEC_1D(self.GIn[it,:,iAnt:iAnt+1,iDir,xPol,yPol],
                                                      self.CentralFreqs,
                                                      Tol=5.e-2,
                                                      Mode=Mode)
                TECMachine.setX0(TEC0CPhase0[:,iAnt].ravel())
                x=TECMachine.doFit()
                LX.append(x)
            X=X1=np.array(LX).reshape((na,len(Mode))).T
                
        if "CPhase" in Mode:
            TEC,CPhase=X.reshape((len(Mode),na))
        else:
            TEC,=X.reshape((len(Mode),na))
            CPhase=np.zeros((1,na),np.float32)
            
        GThis=X2G(X)
        # if "CPhase" in Mode:
        #     TEC,CPhase=X.reshape((len(Mode),na))
        # else:
        #     TEC,=X.reshape((len(Mode),na))
        #     CPhase=np.zeros((1,na),np.float32)
        # TEC-=TEC[0]
        # CPhase-=CPhase[0]
        # GThis=np.abs(GOut[it,:,:,iDir,0,0]).T*TECToZ(TEC.reshape((-1,1)),CPhase.reshape((-1,1)),self.CentralFreqs.reshape((1,-1)))

        # #####################
        # Test
        # GThis1=X2G(X1)

        # print("=======================")
        # print("=======================")
        # print("=======================")
        # print(GThis0)
        # print(GThis)
        # print(GThis1)
        
        # import pylab
        # pylab.clf()
        # pylab.plot(np.angle(GThis).T,color="black")
        # pylab.plot(np.angle(GThis1).T,color="gray")
        # pylab.plot(np.angle(GThis0).T,color="red")
        # pylab.draw()
        # pylab.show()

        
        return GThis.T,TEC,CPhase
        # ######################################
        
        

    
        # # ###########################
        # TEC0,CPhase0=TEC0CPhase0
        # GThis0=TECToZ(TEC0.reshape((-1,1)),CPhase0.reshape((-1,1)),self.CentralFreqs.reshape((1,-1)))

        # for iAnt in range(na):
        #     print "!!!!!!!!!!!!!!",iAnt,iDir
        #     ga=GOut[it,:,iAnt,iDir,0,0]
        #     ge=GThis[iAnt,:]
        #     ge0=GThis0[iAnt,:]
        #     #if iAnt==0: continue
        #     #f=np.linspace(self.CentralFreqs.min(),self.CentralFreqs.max(),100)
        #     #ztec=TECToZ(TECGrid.ravel()[iTec],CPhase.ravel()[iTec],f)
        #     import pylab
        #     pylab.clf()
        #     pylab.subplot(1,2,1)
        #     pylab.scatter(self.CentralFreqs,np.abs(ga),color="black")
        #     pylab.plot(self.CentralFreqs,np.abs(ge),ls=":",color="black")
        #     pylab.subplot(1,2,2)
        #     pylab.scatter(self.CentralFreqs,np.angle(ga),color="black")
        #     pylab.plot(self.CentralFreqs,np.angle(ge),ls=":",color="black")
        #     pylab.plot(self.CentralFreqs,np.angle(ge0),ls=":",color="red")
        #     #pylab.plot(f,np.angle(ztec),ls=":",color="black")
        #     pylab.ylim(-np.pi,np.pi)
        #     pylab.draw()
        #     pylab.show(False)
        #     pylab.pause(0.1)
        # # ###############################
        

       

    def FitThisPolyAmp(self,iDir):
        self._loadDicoStuff()
        for iPol in range(self.NPol):
            xPol,yPol=self.polMap[iPol]
            nt,nch,na,nd,_,_=self.GIn.shape
            GOut=self.GOut # NpShared.GiveArray("%sGOut"%IdSharedMem)
            g=GOut[:,:,:,iDir,xPol,yPol]
            
            AmpMachine=ClassFitAmp.ClassFitAmp(self.GIn[:,:,:,iDir,xPol,yPol],self.CentralFreqs,RemoveMedianAmp=self.RemoveMedianAmp,iDir=iDir)
            gf=AmpMachine.doSmooth()
            #print "Done %i"%iDir
            ga=np.abs(g)
            ga[ga==0]=1
            gf=gf*g/ga
        
            GOut[:,:,:,iDir,xPol,yPol]=gf[:,:,:]
        if self.JonesMode=="Scalar":
            GOut[:,:,:,iDir,1,1]=gf[:,:,:]

    def ClipThisDir(self,iDir,xPol,yPol):
        nt,nch,na,nd,_,_=self.GIn.shape
        GOut=self.GOut#NpShared.GiveArray("%sGOut"%IdSharedMem)
        # g=GOut[:,:,:,iDir,0,0]

        AmpMachine=ClassClip.ClassClip(self.GIn[:,:,:,iDir,xPol,yPol],self.CentralFreqs,RemoveMedianAmp=self.RemoveMedianAmp)
        gf=AmpMachine.doClip()
        GOut[:,:,:,iDir,xPol,yPol]=gf[:,:,:]
        if self.JonesMode=="Scalar":
            GOut[:,:,:,iDir,1,1]=gf[:,:,:]

        
    def FitThisAmp(self,iAnt,iDir):
        self._loadDicoStuff()
        nt,nch,na,nd,_,_=self.GIn.shape

        for iPol in range(self.NPol):
            xPol,yPol=self.polMap[iPol]
            if self.Amp_SmoothType=="Poly":
                for it in range(nt):
                    self.FitThisAmpTimePoly(it,iAnt,iDir,xPol,yPol)
            elif self.Amp_SmoothType=="Gauss":
                self.GaussSmoothAmp(iAnt,iDir,xPol,yPol)

    def EstimateThisTECTime(self,it,iAnt,iDir,xPol,yPol):
        GOut=self.GOut # NpShared.GiveArray("%sGOut"%IdSharedMem)
        g=GOut[it,:,iAnt,iDir,xPol,yPol]
        ga=np.abs(g)
        ga[ga==0]=1
        g0=g/ga
        
        W=np.ones(g0.shape,np.float32)
        W[g==1.]=0
        Z=self.Z
        for iTry in range(5):
            R=(g0.reshape((1,-1))-Z)*W.reshape((1,-1))
            Chi2=np.sum(np.abs(R)**2,axis=1)
            iTec=np.argmin(Chi2)
            rBest=R[iTec]
            if np.max(np.abs(rBest))==0: break
            Sig=np.sum(np.abs(rBest*W))/np.sum(W)
            ind=np.where(np.abs(rBest*W)>5.*Sig)[0]
            if ind.size==0: break
            W[ind]=0

            # gz=TECToZ(TECGrid.ravel()[iTec],CPhase.ravel()[iTec],self.CentralFreqs)
            # import pylab
            # pylab.clf()
            # pylab.subplot(2,1,1)
            # pylab.scatter(self.CentralFreqs,rBest)
            # pylab.scatter(self.CentralFreqs[ind],rBest[ind],color="red")
            # pylab.subplot(2,1,2)
            # pylab.scatter(self.CentralFreqs,rBest)
            # pylab.scatter(self.CentralFreqs[ind],rBest[ind],color="red")
            # pylab.draw()
            # pylab.show()


        
        # # ###########################
        # print iAnt,iDir
        # if iAnt==0: return
        # f=np.linspace(self.CentralFreqs.min(),self.CentralFreqs.max(),100)
        # ztec=TECToZ(TECGrid.ravel()[iTec],CPhase.ravel()[iTec],f)
        # import pylab
        # pylab.clf()
        # pylab.subplot(1,2,1)
        # pylab.scatter(self.CentralFreqs,np.abs(g),color="black")
        # pylab.plot(self.CentralFreqs,np.abs(gz),ls=":",color="black")
        # pylab.plot(self.CentralFreqs,np.abs(gz)-np.abs(g),ls=":",color="red")
        # pylab.subplot(1,2,2)
        # pylab.scatter(self.CentralFreqs,np.angle(g),color="black")
        # pylab.plot(self.CentralFreqs,np.angle(gz),ls=":",color="black")
        # pylab.plot(self.CentralFreqs,np.angle(gz)-np.angle(g),ls=":",color="red")
        # #pylab.plot(f,np.angle(ztec),ls=":",color="black")
        # pylab.ylim(-np.pi,np.pi)
        # pylab.draw()
        # pylab.show(False)
        # pylab.pause(0.1)
        # # ###############################

        t0=self.TECGrid.ravel()[iTec]
        c0=self.CPhase.ravel()[iTec]
    
        gz=np.abs(g)*TECToZ(t0,c0,self.CentralFreqs)
        return gz,t0,c0


        
    def FitThisAmpTimePoly(self,it,iAnt,iDir,xPol,yPol):
        GOut=self.GOut#NpShared.GiveArray("%sGOut"%IdSharedMem)
        g=GOut[it,:,iAnt,iDir,xPol,yPol]
        g0=np.abs(g)
        
        W=np.ones(g0.shape,np.float32)
        W[g0==1.]=0
        # if np.count_nonzero(W)<self.Amp_PolyOrder*3: return
        # stop
        for iTry in range(5):
            if np.max(W)==0: return
            z = np.polyfit(self.CentralFreqs, g0, self.Amp_PolyOrder,w=W)
            # print(it,iAnt,iDir,z)
            p = np.poly1d(z)
            # if iAnt==12:
            #     import pylab
            #     pylab.clf()
            #     pylab.scatter(self.CentralFreqs, g0)
            #     pylab.plot(self.CentralFreqs, p(self.CentralFreqs))
            #     pylab.title(iTry)
            #     pylab.draw()
            #     pylab.show(block=False)
            #     pylab.pause(0.5)
            gz=p(self.CentralFreqs)*g/np.abs(g)
            rBest=(g0-gz)
            if np.max(np.abs(rBest))==0: break
            Sig=np.sum(np.abs(rBest*W))/np.sum(W)
            ind=np.where(np.abs(rBest*W)>5.*Sig)[0]
            if ind.size==0: break
            W[ind]=0

        GOut[it,:,iAnt,iDir,xPol,yPol]=gz
        if self.JonesMode=="Scalar":
            GOut[it,:,iAnt,iDir,1,1]=gz
        

    def GaussSmoothAmp(self,iAnt,iDir,xPol,yPol):
        #print iAnt,iDir
        GOut=self.GOut#NpShared.GiveArray("%sGOut"%IdSharedMem)
        g=GOut[:,:,iAnt,iDir,xPol,yPol]
        g0=np.abs(g)

        
        sg0=scipy.ndimage.filters.gaussian_filter(g0,self.Amp_GaussKernel)

        gz=sg0*g/np.abs(g)
        #print iAnt,iDir,GOut.shape,gz.shape

        GOut[:,:,iAnt,iDir,xPol,yPol]=gz[:,:]
        if self.JonesMode=="Scalar":
            GOut[:,:,iAnt,iDir,1,1]=gz[:,:]
            
        #print np.max(GOut[:,:,iAnt,iDir,0,0]-gz[:,:])

    # def smoothGPR(self):
    #     nt,nch,na,nd,_,_=self.GOut.shape
        
        
        
    def SpacialSmoothTEC(self,iPol,jPol):
        log.print("Do the spacial smoothing...")
        stop
        t=table("/data/tasse/P025+41/L593429_SB132_uv.pre-cal_12A2A9C48t_148MHz.pre-cal.ms/ANTENNA")
        X,Y,Z=t.getcol("POSITION").T
        dx=X.reshape((-1,1))-X.reshape((1,-1))
        dy=Y.reshape((-1,1))-Y.reshape((1,-1))
        dz=Z.reshape((-1,1))-Z.reshape((1,-1))
        D=np.sqrt(dx**2+dy**2+dz**2)
        D0=500.
        WW=np.exp(-D**2/(2.*D0**2))
        WWsum=np.sum(WW,axis=0)
        nt,nch,na,nd,_,_=self.GOut.shape
        
        nt,nd,na = self.TECArray.shape
        for it in range(nt):
            for iDir in range(nd):
                TEC=Tec=self.TECArray[it,iDir]
                TMean=np.dot(WW,Tec.reshape((-1,1))).ravel()
                TMean/=WWsum.ravel()

                # import pylab
                # pylab.clf()
                # pylab.plot(TEC.ravel())
                # pylab.plot(TMean.ravel())
                # pylab.draw()
                # pylab.show(False)
                # pylab.pause(0.5)
                # stop
                

                self.TECArray[it,iDir,:]=TMean[:]

                CPhase=self.CPhaseArray[it,iDir]
                CPMean=np.dot(WW,CPhase.reshape((-1,1))).ravel()
                CPMean/=WWsum.ravel()
                self.CPhaseArray[it,iDir,:]=CPMean[:]

                
                z=np.abs(self.GOut[it,:,:,iDir,0,0]).T*TECToZ(TMean.reshape((-1,1)),
                                                              CPMean.reshape((-1,1)),
                                                              self.CentralFreqs.reshape((1,-1)))
                self.GOut[it,:,:,iDir,0,0]=z.T
                self.GOut[it,:,:,iDir,1,1]=z.T
        
    def Save(self):
        OutFile=self.OutSolsName
        if not ".npz" in OutFile: OutFile+=".npz"

        #self.SpacialSmoothTEC()

        if "TEC" in self.InterpMode:
            # OutFileTEC="%s.TEC_CPhase.npz"%OutFile
            # log.print("  Saving TEC/CPhase solution file as: %s"%OutFileTEC)
            # np.savez(OutFileTEC,
            #          TEC=self.TECArray,
            #          CPhase=self.CPhaseArray)
            self.DicoFile["SolsTEC"]=self.TECArray
            self.DicoFile["SolsCPhase"]=self.CPhaseArray
            

            
        log.print("  Saving interpolated solution file as: %s"%OutFile)
        self.DicoFile["SmoothMode"]=self.InterpMode
        self.DicoFile["SolsOrig"]=copy.deepcopy(self.DicoFile["Sols"])
        self.DicoFile["SolsOrig"]["G"][:]=self.DicoFile["Sols"]["G"][:]
        self.DicoFile["Sols"]["G"][:]=self.GOut[:]
        np.savez(OutFile,**(self.DicoFile))
        
        if "ArraySolsPath" in self.DicoFile.keys():
            log.print("Coming from merged solution file, creating symbolic link")
            for f in self.DicoFile["ArraySolsPath"]:
                if f.decode("ascii")=="": continue
                Dir=os.path.dirname(f).decode("ascii")
                SolName=OutFile.split("/")[-1]
                symsolname="%s/%s"%(Dir,SolName)
                if os.path.islink(symsolname):
                    log.print('Symlink ' + symsolname + ' already exists, recreating')
                    os.unlink(symsolname)
                            
                log.print('Creating ' + symsolname + ' symbolic link')
                os.symlink(os.path.abspath(OutFile),symsolname)
                

            
        # import PlotSolsIm
        # G=self.DicoFile["Sols"]["G"].view(np.recarray)
        # iAnt,iDir=10,0
        # import pylab
        # pylab.clf()
        # A=self.GOut[:,:,iAnt,iDir,0,0]
        # B=G[:,:,iAnt,iDir,0,0]
        # Gs=np.load(OutFile)["Sols"]["G"].view(np.recarray)
        # C=Gs[:,:,iAnt,iDir,0,0]
        # pylab.subplot(1,3,1)
        # pylab.imshow(np.abs(A).T,interpolation="nearest",aspect="auto")
        # pylab.subplot(1,3,2)
        # pylab.imshow(np.abs(B).T,interpolation="nearest",aspect="auto")
        # pylab.subplot(1,3,3)
        # pylab.imshow(np.abs(C).T,interpolation="nearest",aspect="auto")
        # pylab.draw()
        # pylab.show()
        # PlotSolsIm.Plot([self.DicoFile["Sols"].view(np.recarray)])

        NpShared.DelAll("%s"%IdSharedMem)
        #del(self.DicoStuff)

# ############################################        

def test():
    FileName="L401839.killms_f_ap_deep.merged.npz"
    OutFile="TestMerge.Interp.npz"
    CI=ClassInterpol(FileName,OutFile)
    CI.InterpolParallel()
    return CI.Save()

def main(options=None):
    if options==None:
        f = open(SaveName,'rb')
        options = pickle.load(f)
    #FileName="killMS.KAFCA.sols.npz"


    if options.SolsFileIn is None or options.SolsFileOut is None:
        raise RuntimeError("You have to specify In/Out solution file names")
    CI=ClassInterpol(options.SolsFileIn,
                     options.SolsFileOut,
                     InterpMode=options.InterpMode,
                     Amp_PolyOrder=options.Amp_PolyOrder,
                     Amp_GaussKernel=options.Amp_GaussKernel,
                     Amp_SmoothType=options.Amp_SmoothType,
                     NCPU=options.NCPU,CrossMode=options.CrossMode,RemoveMedianAmp=options.RemoveMedianAmp)
    
    #AsyncProcessPool.init(ncpu=NCPU,affinity=0)
    AsyncProcessPool.APP.registerJobHandlers(CI)
    CI.InterpolParallel()
    CI.stopWorkers()
    CI.Save()
    

if __name__=="__main__":
    read_options()
    f = open(SaveName,'rb')
    options = pickle.load(f)
    
    main(options=options)
