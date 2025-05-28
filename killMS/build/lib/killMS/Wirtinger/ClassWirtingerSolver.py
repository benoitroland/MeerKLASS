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
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
#from killMS.Array import NpShared
from DDFacet.Array import shared_dict

from killMS.Predict.PredictGaussPoints_NumExpr5 import ClassPredict
from DDFacet.Other import AsyncProcessPool as APP

from killMS.Data import ClassVisServer
#from Sky import ClassSM
from killMS.Array import ModLinAlg
#import matplotlib.pyplot as pylab

from DDFacet.Other import logger
log=logger.getLogger("ClassWirtingerSolver")
from killMS.Other import ModColor
import os

#from killMS.Other.progressbar import ProgressBar
from DDFacet.Other.progressbar import ProgressBar
            
#from Sky.PredictGaussPoints_NumExpr import ClassPredict
from killMS.Other import ClassTimeIt
from killMS.Other import Counter
from .ClassEvolve import ClassModelEvolution
import time
from itertools import product as ItP
from killMS.Wirtinger import ClassSolverLM
from killMS.Wirtinger import ClassSolverEKF
from killMS.Wirtinger import ClassSolPredictMachine
from SkyModel.Sky import ClassSM
#import killMS.SmoothSols
logger.setSilent("ClassInterpol")
import time



class ClassWirtingerSolver():

    def __init__(self,VS,SM,
                 BeamProps=None,
                 PolMode="IFull",
                 Lambda=1,NIter=20,
                 NCPU=6,
                 SolverType="CohJones",
                 IdSharedMem="",
                 evP_StepStart=0, evP_Step=1,
                 DoPlot=False,
                 DoPBar=True,GD=None,
                 ConfigJacobianAntenna={},
                 TypeRMS="GlobalData",
                 VS_PredictCol=None):
        self.VS_PredictCol=VS_PredictCol
        self.DType=np.complex128
        self.TypeRMS=TypeRMS
        self.IdSharedMem=IdSharedMem
        self.ConfigJacobianAntenna=ConfigJacobianAntenna
        self.Lambda=Lambda
        self.NCPU=NCPU
        self.DoPBar=DoPBar
        self.GD=GD
        self.Q=None
        self.PListKeep=[]
        self.QListKeep=[]
        self.restart_time=0
        
        
        #if self.GD["Solutions"]["SmoothInterpMode"] is not None:
        #    stop
        
        
        # if BeamProps!=None:
        #     rabeam,decbeam=SM.ClusterCat.ra,SM.ClusterCat.dec
        #     Mode,TimeMin=BeamProps
        #     LofarBeam=(Mode,TimeMin,rabeam,decbeam)
        #     VS.SetBeam(LofarBeam)

        MS=VS.MS
        if SM.Type=="Catalog" or SM.Type=="Hybrid":
            SM.Calc_LM(MS.rac,MS.decc)
        self.SM=SM
        self.VS=VS
        self.DoPlot=DoPlot
        if DoPlot==2:
            self.InitPlotGraph()
        self.PolMode=PolMode

        self.SM_Compress=None
        if (self.GD["Compression"]["CompressionMode"] is not None) or self.GD["Compression"]["CompressionDirFile"]:
            log.print(ModColor.Str("Using compression with Mode = %s"%self.GD["Compression"]["CompressionMode"]))
            
            if self.GD["Compression"]["CompressionMode"] and self.GD["Compression"]["CompressionMode"].lower()=="auto":
                ClusterCat=SM.ClusterCat#[1:2]
                #self.SM_Compress=ClassSM.ClassSM(SM)
                
            else:
                ClusterCat=np.load(self.GD["Compression"]["CompressionDirFile"])
                
            ClusterCat=ClusterCat.view(np.recarray)
            SourceCat=np.zeros((ClusterCat.shape[0],),dtype=ClassSM.dtypeSourceList)
            SourceCat=SourceCat.view(np.recarray)
            SourceCat.ra[:]=ClusterCat.ra[:]
            SourceCat.dec[:]=ClusterCat.dec[:]
            SourceCat.RefFreq[:]=100.e6
            SourceCat.I[:]=1.
            SourceCat.Cluster=np.arange(ClusterCat.shape[0])
            np.save("SM_Compress.npy",SourceCat)
            self.SM_Compress=ClassSM.ClassSM("SM_Compress.npy")
            self.SM_Compress.Calc_LM(self.VS.MS.rac,self.VS.MS.decc)
                        
                
        if self.PolMode=="IDiag":
            npolx=2
            npoly=1
        elif self.PolMode=="Scalar":
            npolx=1
            npoly=1
        elif self.PolMode=="IFull":
            npolx=2
            npoly=2

        self.NJacobBlocks_X,self.NJacobBlocks_Y=npolx,npoly

        self.SolPredictMachine=None
        if self.GD["KAFCA"]["EvolutionSolFile"]!="":
            self.SolPredictMachine=ClassSolPredictMachine.ClassSolPredictMachine(GD)
            
        
        self.G=None
        self.NIter=NIter
        #self.SolsList=[]
        self.iCurrentSol=0
        self.SolverType=SolverType
        self.rms=None
        self.rmsFromData=None
        self.SM.ApparentSumI=None
        # if SolverType=="KAFCA":
        #     log.print( ModColor.Str("niter=%i"%self.NIter))
        #     #self.NIter=1
        self.EvolvePStepStart,EvolvePStep=evP_StepStart,evP_Step
        self.CounterEvolveP=Counter.Counter(EvolvePStep)
        self.ThisStep=0
        self.rmsFromExt=None


        
    # def AppendEmptySol(self):
    #     #### Solutions
    #     # self.NSols=self.VS.TimesVisMin.size-1
    #     na=self.VS.MS.na
    #     nd=self.SM.NDir
    #     Sol=np.zeros((1,),dtype=[("t0",np.float64),("t1",np.float64),("G",np.complex64,(na,nd,2,2))])
    #     self.SolsList.append(Sol.view(np.recarray))

    def GiveSols(self,SaveStats=False,OnlyLast=False):
        ind=np.where(self.SolsArray_done==1)[0]
        self.SolsArray_Full.t0[0:ind.size]=self.SolsArray_t0[0:ind.size]
        self.SolsArray_Full.t1[0:ind.size]=self.SolsArray_t1[0:ind.size]
        self.SolsArray_Full.Stats[0:ind.size]=self.SolsArray_Stats[0:ind.size]
        if self.PolMode=="Scalar":
            self.SolsArray_Full.G[0:ind.size,:,:,:,0,0]=self.SolsArray_G[0:ind.size,:,:,:,0,0]
            self.SolsArray_Full.G[0:ind.size,:,:,:,1,1]=self.SolsArray_G[0:ind.size,:,:,:,0,0]
        elif self.PolMode=="IDiag":
            self.SolsArray_Full.G[0:ind.size,:,:,:,0,0]=self.SolsArray_G[0:ind.size,:,:,:,0,0]
            self.SolsArray_Full.G[0:ind.size,:,:,:,1,1]=self.SolsArray_G[0:ind.size,:,:,:,1,0]
        else:                
            self.SolsArray_Full.G[0:ind.size]=self.SolsArray_G[0:ind.size]
            
        Mask=self.SolsArray_HasSolved[...,np.newaxis,np.newaxis,np.newaxis]
        self.SolsArray_Full.HasSolved[0:ind.size,:,:]=self.SolsArray_HasSolved[0:ind.size,:,:]
        
        # if self.GD["Solutions"]["ApplyToDir"]==-2:# and self.GD["Solutions"]["ZeroUnsolved"]:
        #     HasSolved=self.SolsArray_HasSolved
        #     Mask=self.SolsArray_HasSolved[...,np.newaxis,np.newaxis,np.newaxis]
        #     self.SolsArray_Full.G[:]=self.SolsArray_Full.G[:]*Mask
            
        

        # if SaveStats:
        #     ListStd=[l for l in self.ListStd if len(l)>0]
        #     Std=np.array(ListStd)
        #     ListMax=[l for l in self.ListMax if len(l)>0]
        #     Max=np.array(ListMax)
            
        #     ListKapa=[l for l in self.ListKeepKapa if len(l)>0]
        #     Kapa=np.array(ListKapa)
        #     nf,na,nt=Std.shape
        #     NoiseInfo=np.zeros((nf,na,nt,3))
        #     NoiseInfo[:,:,:,0]=Std[:,:,:]
        #     NoiseInfo[:,:,:,1]=np.abs(Max[:,:,:])
        #     NoiseInfo[:,:,:,2]=Kapa[:,:,:]
            
        #     StatFile="NoiseInfo.npy"
        #     log.print( "Saving statistics in %s"%StatFile)
        #     np.save(StatFile,NoiseInfo)

        Sols=self.SolsArray_Full[0:ind.size].copy()
        
            
        if Sols.size==0:
            na=self.VS.MS.na
            nd=self.SM.NDir
            nChan=self.VS.NChanJones
            Sols=np.zeros((1,),dtype=[("t0",np.float64),
                                      ("t1",np.float64),
                                      ("G",np.complex64,(nChan,na,nd,2,2)),
                                      ("Stats",np.float32,(nChan,na,4))])
            Sols=Sols.view(np.recarray)
            Sols.t0[0]=0
            Sols.t1[0]=self.VS.MS.F_times_all[-1]
            Sols.G[0,:,:,:,0,0]=1
            Sols.G[0,:,:,:,1,1]=1
            

        Sols=Sols.view(np.recarray)
        Sols.t1[-1]+=1
        Sols.t0[0]-=1

        return Sols

    def InitSol(self,G=None,TestMode=True):
        na=self.VS.MS.na
        nd=self.SM.NDir
        nChan=self.VS.NChanJones

        if type(G)==type(None):
            if self.PolMode=="Scalar":
                G=np.ones((nChan,na,nd,1,1),self.DType)
            elif self.PolMode=="IDiag":
                G=np.ones((nChan,na,nd,2,1),self.DType)
            else:
                G=np.zeros((nChan,na,nd,2,2),self.DType)
                G[:,:,:,0,0]=1
                G[:,:,:,1,1]=1
            self.HasFirstGuessed=False

        else:
            self.HasFirstGuessed=True
        self.G=G
        self.HasSolvedArray=np.zeros((nChan,na),bool)
        #self.G*=0.001
        _,_,_,npolx,npoly=self.G.shape


        ## # print "int!!!!!!!!!!"
        self.G+=np.random.randn(*self.G.shape)*1.+1j*np.random.randn(*self.G.shape)*1.#sigP
        
        NSols=np.max([1,int(1.5*round(self.VS.MS.DTh/(self.VS.TVisSizeMin/60.)))])
        #print "Nsols",NSols,self.VS.MS.DTh,self.VS.TVisSizeMin/60.
        

        SolsArray_t0=np.zeros((NSols,),dtype=np.float64)
        SolsArray_t1=np.zeros((NSols,),dtype=np.float64)
        SolsArray_tm=np.zeros((NSols,),dtype=np.float64)
        SolsArray_done=np.zeros((NSols,),dtype=np.bool8)
        SolsArray_HasSolved=np.zeros((NSols,nChan,na),dtype=np.bool8)
        SolsArray_G=np.zeros((NSols,nChan,na,nd,npolx,npoly),dtype=np.complex64)
        self.SolsArray_Stats=np.zeros((NSols,nChan,na,4),dtype=np.float32)

        #self.Power0=np.zeros((nChan,na),np.float32)

        
        # self.SolsArray_t0=NpShared.ToShared("%sSolsArray_t0"%self.IdSharedMem,self.SolsArray_t0)
        # self.SolsArray_t1=NpShared.ToShared("%sSolsArray_t1"%self.IdSharedMem,self.SolsArray_t1)
        # self.SolsArray_tm=NpShared.ToShared("%sSolsArray_tm"%self.IdSharedMem,self.SolsArray_tm)
        # self.SolsArray_done=NpShared.ToShared("%sSolsArray_done"%self.IdSharedMem,self.SolsArray_done)
        # self.SolsArray_G=NpShared.ToShared("%sSolsArray_G"%self.IdSharedMem,self.SolsArray_G)
        
        self.DicoSolsArray=shared_dict.create("DicoSolsArray")
        self.DicoSolsArray["t0"]=SolsArray_t0
        self.DicoSolsArray["t1"]=SolsArray_t1
        self.DicoSolsArray["tm"]=SolsArray_tm
        self.DicoSolsArray["done"]=SolsArray_done
        self.DicoSolsArray["G"]=SolsArray_G
        self.DicoSolsArray["HasSolved"]=SolsArray_HasSolved
        self.DicoSolsArray["PredictIsInit"]=np.zeros((NSols,na),np.int16)

        self.SolsArray_t0=self.DicoSolsArray["t0"]
        self.SolsArray_t1=self.DicoSolsArray["t1"]
        self.SolsArray_tm=self.DicoSolsArray["tm"]
        self.SolsArray_done=self.DicoSolsArray["done"]
        self.SolsArray_G=self.DicoSolsArray["G"]
        self.SolsArray_HasSolved=self.DicoSolsArray["HasSolved"]
        
        self.PredictIsInit=self.DicoSolsArray["PredictIsInit"]
        
        self.SolsArray_Full=np.zeros((NSols,),
                                     dtype=[("t0",np.float64),
                                            ("t1",np.float64),
                                            ("G",np.complex64,(nChan,na,nd,2,2)),
                                            ("Stats",np.float32,(nChan,na,4)),
                                            ("HasSolved",bool,(nChan,na))])
        
        
        self.SolsArray_Full=self.SolsArray_Full.view(np.recarray)

        self.DicoKapa={}
        self.DicoKeepKapa={}
        self.DicoStd={}
        self.DicoMax={}
        for (iAnt,iChanSol) in ItP(range(na),range(nChan)):
            self.DicoKapa[(iAnt,iChanSol)]=[]
            self.DicoKeepKapa[(iAnt,iChanSol)]=[]
            self.DicoStd[(iAnt,iChanSol)]=[]
            self.DicoMax[(iAnt,iChanSol)]=[]


        self.DicoCurrentGains=shared_dict.create("DicoCurrentGains")
        self.DicoCurrentGains["G"]=self.G.copy()
        self.DicoCurrentGains["G1"]=self.G.copy()
        self.DicoCurrentGains["G0Iter"]=self.G.copy()
        self.G=self.DicoCurrentGains["G"]
        self.G1=self.DicoCurrentGains["G1"]
        self.G0Iter=self.DicoCurrentGains["G0Iter"]
        self.DicoCurrentGains["FreqsDone"]=np.zeros((na,nChan),np.int16)
        
        # self.G=NpShared.ToShared("%sSharedGains"%self.IdSharedMem,self.G)
        # self.G0Iter=NpShared.ToShared("%sSharedGains0Iter"%self.IdSharedMem,self.G.copy())
        
        
        #self.InitCovariance()

    def GiveDicoSols(self,Sols=None,OnlyLast=False):
        if Sols is None:
            Sols=self.GiveSols(OnlyLast=OnlyLast)
        SolsSave=Sols
        ClusterCat=self.SM.ClusterCat
        nch,na,nd,_,_=self.G.shape
        
        MeanDataFreqPerSolBand=(self.VS.SolsFreqDomains[:,0]+self.VS.SolsFreqDomains[:,1])/2.
        
        if self.SM.Type=="Image" or self.SM.Type=="Catalog" or self.SM.Type=="Hybrid":
            nt,nch,na,nd,_,_=Sols.G.shape
            nd=self.SM.NDirsOrig
            SolsAll=np.zeros((nt,),dtype=[("t0",np.float64),("t1",np.float64),("G",np.complex64,(nch,na,nd,2,2)),("Stats",np.float32,(nch,na,4))])
            SolsAll=SolsAll.view(np.recarray)
            SolsAll.G[:,:,:,:,0,0]=1
            SolsAll.G[:,:,:,:,1,1]=1
            SolsAll.G[:,:,:,self.SM.MapClusterCatOrigToCut,:,:]=Sols.G[:,:,:,:,:,:]
            SolsAll.t0=Sols.t0
            SolsAll.t1=Sols.t1
            SolsAll.Stats=Sols.Stats
            SolsSave=SolsAll
            ClusterCat=self.SM.ClusterCatOrig
            try:
                if self.SM.Type=="Hybrid":
                    MeanDataFreqPerSolBand=self.LSM[0].MeanDataFreqPerSolBand
                else:
                    MeanDataFreqPerSolBand=self.SM.MeanDataFreqPerSolBand
                    
            except:
                MeanDataFreqPerSolBand=(self.VS.SolsFreqDomains[:,0]+self.VS.SolsFreqDomains[:,1])/2.
            
        StationNames=np.array(self.VS.MS.StationNames)

        nPreApplyTimes=self.VS.PreApplyTimesT0.size
        PreApplyTimes=np.zeros((nPreApplyTimes,),dtype=[("t0",np.float64),("t1",np.float64)]).view(np.recarray)
        PreApplyTimes.t0=self.VS.PreApplyTimesT0
        PreApplyTimes.t1=self.VS.PreApplyTimesT1

        nBeamTimes=self.VS.BeamTimesT0.size
        BeamTimes=np.zeros((nBeamTimes,),dtype=[("t0",np.float64),("t1",np.float64)]).view(np.recarray)
        BeamTimes.t0=self.VS.BeamTimesT0
        BeamTimes.t1=self.VS.BeamTimesT1
        
        D={"MSName":os.path.abspath(self.VS.MS.MSName),
           "MSNameTime0":self.VS.MS.Time0,
           "Sols":SolsSave,
           "StationNames":StationNames,
           "SkyModel":ClusterCat,
           "ClusterCat":ClusterCat,
           #"SourceCatSub":SourceCatSub,
           "ModelName":self.GD["SkyModel"]["SkyModel"],
           "FreqDomains":self.VS.SolsFreqDomains,
           "PreApplyTimes":PreApplyTimes,
           "BeamTimes":BeamTimes,
           "MeanDataFreqPerSolBand":MeanDataFreqPerSolBand}
        return D
        
    # def InitSmoother(self):
    #     if self.GD["Solutions"]["SmoothInterpMode"] is not None:
    #         log.print(ModColor.Str("Will smooth solutions using %s mode"%self.GD["Solutions"]["SmoothInterpMode"],col="green"))
    #         D=self.GiveDicoSols()
    #         self.CI=killMS.SmoothSols.ClassInterpol(InSolsName=D,
    #                                                 OutSolsName=None,
    #                                                 InterpMode=self.GD["Solutions"]["SmoothInterpMode"],
    #                                                 DoPBar=False,
    #                                                 NCPU=self.NCPU)
    #         self.CI.startWorkers()
            
    # def setSmoother(self,CI):
    #     self.CI=CI
    #     D=self.GiveDicoSols()
    #     self.CI.setInSols(D)
    #     #self.CI.startWorkers()
        
    # def SmoothSols(self):
    #     if not self.GD["Smoother"]["SmoothEnable"]: return
    #     na=self.VS.MS.na
    #     nd=self.SM.NDir
    #     nChan=self.VS.NChanJones
    #     Sols=np.zeros((1,),dtype=[("t0",np.float64),
    #                               ("t1",np.float64),
    #                               ("G",np.complex64,(nChan,na,nd,2,2)),
    #                               ("Stats",np.float32,(nChan,na,4))])
    #     Sols=Sols.view(np.recarray)
    #     Sols.t0[0]=0
    #     Sols.t1[0]=self.VS.MS.F_times_all[-1]

    #     if self.PolMode=="Scalar":
    #         Sols.G[0,:,:,:,0,0]=self.G[:,:,:,0,0]
    #         Sols.G[0,:,:,:,1,1]=self.G[:,:,:,0,0]
    #     elif self.PolMode=="IDiag":
    #         Sols.G[0,:,:,:,0,0]=self.G[:,:,:,0,0]
    #         Sols.G[0,:,:,:,1,1]=self.G[:,:,:,1,0]
    #     else:                
    #         Sols.G[0,:,:,:,:,:]=self.G[:,:,:,:,:]
            
    #     D=self.GiveDicoSols(Sols=Sols)
    #     self.CI.updateSols(D)
    #     self.CI.InterpolParallel()

    #     gg=self.G[:,:,:,0,0].copy()
    #     if self.PolMode=="Scalar":
    #         GO=self.CI.GOut[0,:,:,:,0,0].copy()
    #         self.G[:,:,:,0,0]=GO[:,:,self.SM.MapClusterCatOrigToCut]
    #     else:                
    #         stop


        
    def InitCovariance(self,FromG=False,sigP=0.1,sigQ=0.01):
        if self.SolverType!="KAFCA": return
        if self.Q!=None: return
        
        na=self.VS.MS.na
        nd=self.SM.NDir
        nChan=self.VS.NChanJones

        
        _,_,_,npol,_=self.G.shape
        


        npolx,npoly=self.NJacobBlocks_X,self.NJacobBlocks_Y

        if FromG==False:
            P=(sigP**2)*np.array([np.diag(np.ones((nd*npolx*npoly,),self.DType)) for iAnt in range(na)])
            Q=(sigQ**2)*np.array([np.diag(np.ones((nd*npolx*npoly,),self.DType)) for iAnt in range(na)])
        else:

            P=(sigP**2)*np.array([np.max(np.abs(self.G[iAnt]))**2*np.diag(np.ones((nd*npolx*npoly),self.DType)) for iAnt in range(na)])
            Q=(sigQ**2)*np.array([np.max(np.abs(self.G[iAnt]))**2*np.diag(np.ones((nd*npolx*npoly),self.DType)) for iAnt in range(na)])

        if self.SM.ApparentSumI is None:
            self.InitMeanBeam()


        #self.VS.giveDataSizeAntenna()
        QList=[]
        PList=[]
        for iChanSol in range(self.VS.NChanJones):
            ra=self.SM.ClusterCat.ra
            dec=self.SM.ClusterCat.dec

            ns=ra.size
            
            d=np.sqrt((ra.reshape((ns,1))-ra.reshape((1,ns)))**2+(dec.reshape((ns,1))-dec.reshape((1,ns)))**2)
            d0=1.*np.pi/180
            QQ=(1./(1.+d/d0))**2
            Qa=np.zeros((nd,npolx,npoly,nd,npolx,npoly),self.DType)
            for ipol in range(npolx):
                for jpol in range(npoly):
                    Qa[:,ipol,jpol,:,ipol,jpol]=QQ[:,:]

            #Qa=np.zeros((nd,npolx,npoly,nd,npolx,npoly),self.DType)
            F=self.SM.ClusterCat.SumI.copy()
            F/=F.max()

            #stop
            #self.SM.ApparentSumI=np.zeros((nd,),np.float32)
            Qa.fill(0)
            ApFluxes=self.NormFluxes*self.AbsMeanBeamAnt**2
            
            for idir in range(nd):
                #Qa[idir,:,:,idir,:,:]*=(self.SM.ApparentSumI[idir])**2
                #Qa[idir,:,:,idir,:,:]*=ApFluxes[idir]**2
                #Qa[idir,:,:,idir,:,:]=1
                Qa[idir,:,:,idir,:,:]=ApFluxes[idir]**2


            # for idir in range(nd):
            #     for jdir in range(nd):
            #         Qa[idir,:,:,jdir,:,:]=np.sqrt(ApFluxes[idir]*ApFluxes[jdir])*QQ[idir,jdir]

            # import pylab
            # pylab.clf()
            # pylab.imshow(QQ,interpolation="nearest")
            # pylab.draw()
            # pylab.show()

            
            Qa=Qa.reshape((nd*npolx*npoly,nd*npolx*npoly))
            #print np.diag(Qa)
            #Q=(sigQ**2)*np.array([np.max(np.abs(self.G[iChanSol,iAnt]))**2*(Qa*(self.VS.fracNVisPerAnt[iAnt]**4))**(self.GD["KAFCA"]["PowerSmooth"]) for iAnt in range(na)])
            Q=(sigQ**2)*np.array([np.max(np.abs(self.G[iChanSol,iAnt]))**2*Qa for iAnt in range(na)])
            #Q=(sigQ**2)*np.array([np.max(np.abs(self.G[iChanSol,iAnt]))**2*(Qa*(self.VS.Compactness[iAnt]**2*self.VS.fracNVisPerAnt[iAnt]**4))**(self.GD["KAFCA"]["PowerSmooth"]) for iAnt in range(na)])
            #print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            #Q=(sigQ**2)*np.array([np.max(np.abs(self.G[iChanSol,iAnt]))**2*Qa for iAnt in range(na)])
            #print Q[0]

            # ###############################

            # Qa=np.zeros((na,nd,npolx,npoly,nd,npolx,npoly),self.DType)
            # Qa.fill(0)
            # ApFluxes=self.NormFluxes.reshape((-1,1))*self.AbsMeanBeam**2
            # for iAnt in range(na):
            #     for idir in range(nd):
            #         Qa[iAnt,idir,:,:,idir,:,:]=ApFluxes[idir,iAnt]**0.5#**2

            # Qa=Qa.reshape((na,nd*npolx*npoly,nd*npolx*npoly))
            # Q=(sigQ**2)*np.array([np.max(np.abs(self.G[iChanSol,iAnt]))**2*Qa[iAnt] for iAnt in range(na)])

            QList.append(Q)
            PList.append(P)
        Q=np.array(QList)
        P=np.array(PList)

        
        evP=np.zeros_like(P)
        Q_Init=Q.copy()
        self.DicoCov=shared_dict.create("DicoCov")
        self.DicoCov["P"]=P
        self.DicoCov["Q"]=Q
        self.DicoCov["evP"]=evP
        self.DicoCov["Q_Init"]=Q_Init
        self.P=self.DicoCov["P"]
        self.Q=self.DicoCov["Q"]
        self.Q_Init=self.DicoCov["Q_Init"]
        self.evP=self.DicoCov["evP"]
        
        # self.P=NpShared.ToShared("%sSharedCovariance"%self.IdSharedMem,self.P)
        # self.Q=NpShared.ToShared("%sSharedCovariance_Q"%self.IdSharedMem,Q)
        # self.evP=NpShared.ToShared("%sSharedEvolveCovariance"%self.IdSharedMem,self.evP)
        
        nbuff=10

    def InitMeanBeam(self):

        self.NormFluxes=self.SM.ClusterCat.SumI.copy()
        # print np.sort(self.NormFluxes)
        # FCut=5.
        # self.NormFluxes[self.NormFluxes>FCut]=FCut
        #stop
        isATeam=[("ATeam" in name.decode("ascii")) for name in self.SM.ClusterCat.Name]
        isNotATeam=[("ATeam" not in name.decode("ascii")) for name in self.SM.ClusterCat.Name]
        
        MaxF=self.NormFluxes[isNotATeam].max()
        iBrightest=np.where(self.NormFluxes==MaxF)[0]
        
        self.NormFluxes/=MaxF
        self.NormFluxes[isATeam]=1.
        
        if self.GD["Beam"]["BeamModel"] is None:
            self.SM.ApparentSumI=self.NormFluxes
            self.SM.AbsMeanBeamAnt=np.ones_like(self.SM.ApparentSumI)
            self.AbsMeanBeamAnt=self.SM.AbsMeanBeamAnt
            self.AbsMeanBeam=np.ones((self.SM.ApparentSumI.shape[0],self.VS.MS.na),float)
            self.SM.AbsMeanBeam=self.AbsMeanBeam
        else:
            nd=self.SM.ClusterCat.SumI.size
            self.SM.ApparentSumI=np.zeros((nd,),np.float32)
            from killMS.Data import ClassBeam
            log.print( "Calculate mean beam for covariance estimate... ")
            BeamMachine=ClassBeam.ClassBeam(self.VS.MSName,self.GD,self.SM)
            AbsMeanBeam=BeamMachine.GiveMeanBeam()
            
            AbsMeanBeamAnt=np.mean(AbsMeanBeam[:,:,0,0,0],axis=1)
            
            self.AbsMeanBeam=AbsMeanBeam[:,:,0,0,0]

            self.AbsMeanBeamAnt=AbsMeanBeamAnt
            self.SM.ApparentSumI=(AbsMeanBeamAnt)*self.NormFluxes
            self.SM.AbsMeanBeamAnt=AbsMeanBeamAnt


        for iDir,okATeam in enumerate(isATeam):
            if not okATeam: continue
            self.AbsMeanBeam[iDir,:]=self.AbsMeanBeam[iBrightest,:]
            self.AbsMeanBeamAnt[iDir]=self.AbsMeanBeamAnt[iBrightest]
            self.SM.ApparentSumI=self.AbsMeanBeamAnt*self.NormFluxes
            self.SM.AbsMeanBeamAnt=self.AbsMeanBeamAnt

        # pylab.clf()
        # pylab.scatter(self.SM.ClusterCat.l,self.SM.ClusterCat.m,c=self.SM.ApparentSumI)
        # pylab.draw()
        # pylab.show(False)
        # pylab.pause(0.1)
        # stop

        self.InitReg()

    def InitReg(self):
        if self.SolverType=="KAFCA": return
        if self.GD["CohJones"]["LambdaTk"]==0: return
        NDir=self.SM.NDir
        X0=np.ones((NDir,self.NJacobBlocks_X,self.NJacobBlocks_Y),dtype=np.float32)
        L=np.ones((NDir,self.NJacobBlocks_X,self.NJacobBlocks_Y),dtype=np.float32)
        if self.PolMode=="IFull":
            X0[:,0,1]=0
            X0[:,1,0]=0
                
        SumI=self.SM.ClusterCat.SumI.copy()
        SumIApp=SumI*self.SM.AbsMeanBeamAnt**2

        MaxFlux=1.
        #indFree=np.where(SumIApp>MaxFlux)[0]
        #SumIApp[indFree]=MaxFlux

        #SumIAppNorm=SumIApp#/MaxFlux
        #Linv=L/SumIAppNorm.reshape((NDir,1,1))

        AbsMeanBeamAntsq=self.SM.AbsMeanBeamAnt**2
        Linv=L/AbsMeanBeamAntsq.reshape((NDir,1,1))


        Linv=Linv**2

        log.print( "Using Tikhonov regularisation [LambdaTk = %.2f]"%self.GD["CohJones"]["LambdaTk"])
        #log.print( "  there are %i free directions"%indFree.size)
        log.print( "  minimum inverse L-matrix is %.3f"%Linv.min())
        log.print( "  maximum inverse L-matrix is %.3f"%Linv.max())
        # for iDir in range(NDir):
        #     log.print( "  #%i : [%7.3fJy x %7.7f] %7.3f Jy -> %7.7f "%(iDir,SumI[iDir],self.SM.AbsMeanBeamAnt[iDir],SumIApp[iDir],Linv.flat[iDir]))

        self.DicoTikhonov=shared_dict.create("DicoTikhonov")
        self.DicoTikhonov["Linv"]=Linv
        self.DicoTikhonov["X0"]=X0
        #NpShared.ToShared("%sLinv"%self.IdSharedMem,Linv)
        #NpShared.ToShared("%sX0"%self.IdSharedMem,X0)



        
    def setNextData(self):

        
            

        
        DATA=self.VS.GiveNextVis()
        if self.VS_PredictCol is not None:
            self.DataPredictCol=self.VS_PredictCol.GiveNextVis()

        # print(DATA)

        if DATA=="EndOfObservation":
            log.print( ModColor.Str("Reached end of data"))
            return "EndOfObservation"
        elif DATA=="EndChunk":
            log.print( ModColor.Str("Reached end of data chunk"))
            return "EndChunk"
        elif DATA=="AllFlaggedThisTime":
            #log.print( ModColor.Str("AllFlaggedThisTime"))
            self.AppendGToSolArray()
            self.iCurrentSol+=1
            return "AllFlaggedThisTime"

        
        ## simul
        #d=self.DATA["data"]
        #self.DATA["data"]+=(self.rms/np.sqrt(2.))*(np.random.randn(*d.shape)+1j*np.random.randn(*d.shape))
        self.DATA=DATA

        self.rms=-1
        if (self.TypeRMS=="Resid")&(self.rmsFromData!=None):
            self.rms=self.rmsFromData
            #log.print(" rmsFromDataJacobAnt: %s"%self.rms)
        elif self.rmsFromExt!=None:
            self.rms=self.rmsFromExt
            #log.print(" rmsFromExt: %s"%self.rms)
        elif (self.TypeRMS=="GlobalData"):
            nrow,nch,_=DATA["flags"].shape
            if self.VS.MS.NPolOrig==4:
                Dpol=DATA["data"][:,:,1:3]
                Fpol=DATA["flags"][:,:,1:3]
                w=DATA["W"].reshape((nrow,nch,1))*np.ones((1,1,2))
            else:
                Dpol=DATA["data"][:,:,0:1]
                Fpol=DATA["flags"][:,:,0:1]
                w=DATA["W"].reshape((nrow,nch,1))
                
            self.rms=np.sqrt(np.sum((w[Fpol==0]*np.absolute(Dpol[Fpol==0]))**2.0)/np.sum(w[Fpol==0]**2.0))/np.sqrt(2.)
            # print
            # log.print(" rmsFromGlobalData: %s"%self.rms)
            # print DATA["data"].shape
            # print
        else:
            stop

        
        #print("rms=",self.rms)

        return True

    def SetRmsFromExt(self,rms):
        self.rmsFromExt=rms

    def InitPlotGraph(self):
        from Plot import Graph
        log.print("Initialising plots ..." )
        import pylab
        #pylab.ion()
        self.Graph=Graph.ClassMplWidget(self.VS.MS.na)
        
        for iAnt in range(self.VS.MS.na):
            self.Graph.subplot(iAnt)
            self.Graph.imshow(np.zeros((10,10),dtype=np.float32),interpolation="nearest",aspect="auto",origin='lower',vmin=0.,vmax=2.)#,extent=(-3,3,-3,3))
            self.Graph.text(0,0,self.VS.MS.StationNames[iAnt])
            self.Graph.draw()

        pylab.draw()
        pylab.show(False)


    def AppendGToSolArray(self):
        t0,t1=self.VS.CurrentVisTimes_MS_Sec
        self.SolsArray_t0[self.iCurrentSol]=t0
        self.SolsArray_t1[self.iCurrentSol]=t1
        tm=(t0+t1)/2.
        self.SolsArray_tm[self.iCurrentSol]=tm
        self.SolsArray_done[self.iCurrentSol]=1
        self.SolsArray_G[self.iCurrentSol][:]=self.G[:]
        self.SolsArray_HasSolved[self.iCurrentSol][:]=self.HasSolvedArray[:]
        
        # if self.SolverType=="KAFCA":
        #     self.PListKeep.append(self.P.copy())
        #     self.QListKeep.append(self.Q.copy())
        # NDone=np.count_nonzero(self.SolsArray_done)
        # print(NDone)
        # print(NDone)
        # print(NDone)
        # if NDone>=867:
        #     FileName="CurrentSols.%i.npy"%NDone
        #     log.print( "Save Solutions in file: %s"%FileName)
        #     log.print( "Save Solutions in file: %s"%FileName)
        #     log.print( "Save Solutions in file: %s"%FileName)
        #     Sols=self.GiveSols()
        #     np.save(FileName,Sols.G.copy())
        #     if self.SolverType=="KAFCA":
        #         # np.save("P.%s.npy"%self.GD["Solutions"]['OutSolsName'],np.array(self.PListKeep))
        #         FName="P.%s.%i.npz"%(self.GD["Solutions"]['OutSolsName'],NDone)
        #         log.print( "Save PQ in file: %s"%FName)
        #         log.print( "Save PQ in file: %s"%FName)
        #         log.print( "Save PQ in file: %s"%FName)
        #         np.savez(FName,
        #                  ListP=np.array(self.PListKeep),
        #                  ListQ=np.array(self.QListKeep))
        
        # if self.SolverType=="KAFCA":
        #     self.PListKeep.append(self.P.copy())
        #     self.QListKeep.append(self.Q.copy())
        #     np.save("P.%s.npy"%self.GD["Solutions"]['OutSolsName'],np.array(self.PListKeep))
        #     np.savez("P.%s.npz"%self.GD["Solutions"]['OutSolsName'],
        #              ListP=np.array(self.PListKeep),
        #              ListQ=np.array(self.QListKeep))




    #======================================

    def print_worker(self):
        for iWorker in range(APP.APP.ncpu):
            APP.APP.runJob("initPM_worker:%i"%(iWorker,),
                           self._print_worker,
                           args=(0,))
            
        workers_res=APP.APP.awaitJobResults("initPM_worker*")#, progress=1)

    def _print_worker(self,a):
        print(self.SM.MapClusterCatOrigToCut)
        print(self.PM)
        
        import time
        time.sleep(1)

    def setSM(self,SM):
        self.SM=SM
        
    def InitPM(self):
        
        #self.SM._model_dict=shared_dict.attach(self.SM._model_dict_path)
        
        x=np.linspace(0.,15,100000)
        Exp=np.float32(np.exp(-x))
        LExp=[Exp,x[1]-x[0]]

        Sinc=np.zeros(x.shape,np.float32)
        Sinc[0]=1.
        Sinc[1::]=np.sin(x[1::])/(x[1::])
        LSinc=[Sinc,x[1]-x[0]]

        self.PM=ClassPredict(Precision="S",DoSmearing=self.GD["SkyModel"]["Decorrelation"],
                             IdMemShared=self.IdSharedMem,
                             LExp=LExp,LSinc=LSinc,
                             BeamAtFacet=(self.GD["Beam"]["BeamAt"].lower() == "facet"))

        if self.GD["ImageSkyModel"]["BaseImageName"]!="" and self.GD["SkyModel"]["SkyModelCol"] is None and self.GD["FITSSkyModel"]["FITSSkyModel"] is None:
            self.PM.InitGM(self.SM)

        self.PM_Compress=None
        if self.SM_Compress:
            self.PM_Compress=ClassPredict(Precision="S",
                                          DoSmearing=self.GD["SkyModel"]["Decorrelation"],
                                          IdMemShared=self.IdSharedMem,
                                          LExp=LExp,LSinc=LSinc)
        time.sleep(1)


