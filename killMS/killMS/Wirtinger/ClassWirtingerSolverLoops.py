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
import killMS.SmoothSols
logger.setSilent("ClassInterpol")
import time
from killMS.Other import ClassFitMultiModel

from killMS.Wirtinger.ClassWirtingerSolver import ClassWirtingerSolver

    
class ClassWirtingerSolverLoops(ClassWirtingerSolver):

    def __init__(self,*args,**kwargs):
        ClassWirtingerSolver.__init__(self,*args,**kwargs)
        self.restart_time = time.time()
        
    # #################################
    # ###        Parallel           ###
    # #################################
    
    def doNextTimeSolve_Parallel(self,OnlyOne=False,SkipMode=False,Parallel=True):
        
        Parallel=True
        
        # Parallel=False
        # SkipMode=True
        
        ListAntSolve=[i for i in range(self.VS.MS.na) if not(i in self.VS.FlagAntNumber)]
        #print(ListAntSolve)

        NCPU=self.NCPU

        import time

        #T=ClassTimeIt.ClassTimeIt()
        #T.disable()

        
        ##############################

        T0,T1=self.VS.TimeMemChunkRange_sec[0],self.VS.TimeMemChunkRange_sec[1]
        DT=(T1-T0)
        dt=self.VS.TVisSizeMin*60.
        dt=np.min([dt,DT])
        nt=int(np.ceil(DT/float(dt)))
        nt=self.VS.NSlicePerChunk[self.VS.iChunkThis]
        # if DT/float(dt)-nt>1.:
        #     nt+=1
        #nt=np.max([1,nt])
        
        log.print("DT=%f, dt=%f, nt=%f"%(DT,dt,nt))
        

        Title="Solving [Chunk %i/%i] "%(self.VS.iChunkThis+1,self.VS.NTChunk)

        self.pBAR= ProgressBar(Title=Title)
        #self.pBAR.disable()
        if not(self.DoPBar): self.pBAR.disable()
        
        self.pBAR.render(0,nt)
        NDone=0
        iiCount=0
        ThisG=self.G.copy()
        if self.SolverType=="KAFCA":
            ThisP=self.P.copy()
            ThisQ=self.Q.copy()
            
        while True:
            T=ClassTimeIt.ClassTimeIt("ClassWirtinger DATA[%4.4i]"%NDone)
            T2=ClassTimeIt.ClassTimeIt("ClassWirtingerFull")
            T.disable()
            T2.disable()
            #print("===============================================================")
            T.reinit()
            NDone+=1
            self.pBarProgress=NDone,float(nt)
            Res=self.setNextData()
            T.timeit("read data")
            if Res=="EndChunk":
                break
            if Res=="AllFlaggedThisTime":
                continue
            #print "saving"
            #print "saving"
            #sols=self.GiveSols()
            #np.save("lastSols",sols)
            #print "done"
            if SkipMode:
                continue
                if iiCount<=21:
                    print(iiCount)
                    print(iiCount)
                    iiCount+=1
                    continue

            iiCount+=1
            # print("iiCount",self.VS.iChunkThis,iiCount)
            # # if self.VS.iChunkThis!=4: continue
            # # if iiCount<7: continue
            # print("   doSoleve")
            

            t0,t1=self.VS.CurrentVisTimes_MS_Sec
            tm=(t0+t1)/2.
  
            NJobs=len(ListAntSolve)
            NTotJobs=NJobs*self.NIter

            lold=0
            iResult=0

            T.timeit("stuff")

            if (not(self.HasFirstGuessed))&(self.SolverType=="CohJones"):
                NIter=15
                self.HasFirstGuessed=True
            else:
                NIter=self.NIter


            #print "!!!!!!!!!!!!!!!!!!!!!!!!!"
            #NIter=1


            Gold=self.G.copy()
            DoCalcEvP=np.zeros((self.VS.NChanJones,),bool)            
            DoCalcEvP[:]=False
            if (self.CounterEvolveP())&(self.SolverType=="KAFCA")&(self.iCurrentSol>self.EvolvePStepStart):
                DoCalcEvP[:]=True
            elif (self.SolverType=="KAFCA")&(self.iCurrentSol<=self.EvolvePStepStart):
                DoCalcEvP[:]=True

            T.timeit("before iterloop")
            u,v,w=self.DATA["uvw"].T
            A0=self.DATA["A0"]
            A1=self.DATA["A1"]
            meanW=np.zeros((self.VS.MS.na,),np.float32)
            for iAntMS in ListAntSolve:
                ind=np.where((A0==iAntMS)|(A1==iAntMS))[0]
                if ind.size==0: continue
                meanW[iAntMS]=np.mean(np.abs(w[ind]))
            meanW=meanW[ListAntSolve]
            indOrderW=np.argsort(meanW)[::-1]
            SortedWListAntSolve=(np.array(ListAntSolve)[indOrderW]).tolist()
            # print(SortedWListAntSolve)

            #print indOrderW
            #NpShared.ToShared("%sSharedGainsPrevious"%self.IdSharedMem,self.G.copy())
            #NpShared.ToShared("%sSharedPPrevious"%self.IdSharedMem,self.P.copy())
            Dico_SharedDicoDescriptors={}

            if self.SolPredictMachine is not None:
                t0_ms,t1_ms=self.VS.CurrentVisTimes_MS_Sec
                tm_ms=(t0_ms+t1_ms)/2.
                xPredict=self.SolPredictMachine.GiveClosestSol(tm_ms,
                                                               self.VS.SolsFreqDomains,np.arange(self.VS.MS.na),
                                                               self.SM.ClusterCat.ra,self.SM.ClusterCat.dec)
                
                # nChan,na,nd,2,2
                if self.PolMode=="Scalar":
                    self.G[:,:,:,0,0]=xPredict[:,:,:,0,0]
                elif self.PolMode=="IDiag":
                    self.G[:,:,:,0,0]=xPredict[:,:,:,0,0]
                    self.G[:,:,:,1,0]=xPredict[:,:,:,1,1]
                else:
                    self.G[:]=xPredict[:]


            #self.PredictIsInit.fill(0)
            
            T.timeit("before Init predict")
            # Initialise Predict
            DoEvP=np.zeros((self.VS.NChanJones,),bool)
            LMIter=0
            iChanSol=0
            ThisG[iChanSol,:]=self.G[iChanSol,:]
            if self.SolverType=="KAFCA":
                ThisP[iChanSol,:]=self.P[iChanSol,:]
                ThisQ[iChanSol,:]=self.Q[iChanSol,:]
            self.G0Iter[iChanSol,:]=ThisG[iChanSol,:]
            DoEvP[iChanSol]=False
            DoFullPredict=False



            # ######################################
            if self.GD["Solvers"]["AsyncPredict"]:
                DataDesc_iSliceSolveThis=self.VS.DicoDataDesc.get(self.VS.iSliceSolveThis,None)
                iChunkThisSlice=None
                if DataDesc_iSliceSolveThis is not None:
                    iChunkThisSlice=DataDesc_iSliceSolveThis["iChunk"]
                    
                isFirstSliceOfChunk=self.VS.DicoDataDesc[self.VS.iSliceSolveThis]["isFirstSliceOfChunk"]
                if isFirstSliceOfChunk and (iChunkThisSlice==self.VS.iChunkThis):
                    for iAnt in SortedWListAntSolve:
                        args=(iAnt,
                              self.rms,
                              (self.VS.iSliceSolveThis,self.VS.iChunkThis))
                        APP.APP.runJob("_worker_initJM:Chunk=%i,Slice=%i,Ant=%i"%(self.VS.iChunkThis,self.VS.iSliceSolveThis,iAnt),
                                       self._worker_initJM,
                                       args=args)#,serial=True)
                    
                DataDesc_iSliceSolveNext=self.VS.DicoDataDesc.get(self.VS.iSliceSolveNext,None)
                iChunkNextSlice=None
                if DataDesc_iSliceSolveNext is not None:
                    iChunkNextSlice=DataDesc_iSliceSolveNext["iChunk"]
                if (iChunkNextSlice==self.VS.iChunkThis):
                    for iAnt in SortedWListAntSolve:
                        args=(iAnt,
                              self.rms,
                              (self.VS.iSliceSolveNext,self.VS.iChunkThis))
                        APP.APP.runJob("_worker_initJM:Chunk=%i,Slice=%i,Ant=%i"%(self.VS.iChunkThis,self.VS.iSliceSolveNext,iAnt),
                                       self._worker_initJM,
                                       args=args)#,serial=True)
                
            # ######################################
            
            for LMIter in range(NIter):
                self.HasSolvedArray.fill(0)
                T3=ClassTimeIt.ClassTimeIt("  iter[%i]"%LMIter)
                T3.disable()
                isLastLMIter=(LMIter==(NIter-1))
                # # Reset Data
                # # _,na,_,_,_=self.G.shape
                # g=self.G[iChanSol,:,:,0,0]
                # P=np.mean(np.abs(g)**2,1)
                # if self.iCurrentSol!=0:
                #     fact=self.Power0[iChanSol]/P
                #     self.G[iChanSol]*=fact.reshape((na,1,1,1))
                #     print fact,self.Power0[iChanSol],P
                # else:
                #     self.Power0[iChanSol]=P
                #     print self.Power0[iChanSol]

                #NpShared.DelAll("%sDicoData"%self.IdSharedMem)
                T.reinit()
                self.DicoCurrentGains["FreqsDone"].fill(0)
                T3.timeit("bef")
                for iChanSol in range(self.VS.NChanJones):
                    for iAnt in SortedWListAntSolve:

                        ThisG[iChanSol,:]=self.G[iChanSol,:]

                        if self.SolverType=="KAFCA":
                            ThisP[iChanSol,:]=self.P[iChanSol,:]
                            ThisQ[iChanSol,:]=self.Q[iChanSol,:]

                        #print
                        if LMIter>0:
                            DoCalcEvP[iChanSol]=False
                            DoEvP[iChanSol]=False
                        elif LMIter==0:
                            self.G0Iter[iChanSol,:]=ThisG[iChanSol,:]
                            DoEvP[iChanSol]=False
    
                        if LMIter==(NIter-1):
                            DoEvP[iChanSol]=True
                    
                        DoFullPredict=False
                        #if LMIter==(NIter-1):
                        #    DoFullPredict=True
                        
                        args=(iAnt,iChanSol,
                              DoCalcEvP[iChanSol],tm,
                              self.rms,
                              DoEvP[iChanSol],
                              DoFullPredict,
                              (self.VS.iSliceSolveThis,self.VS.iChunkThis),
                              iiCount,LMIter)
                        APP.APP.runJob("worker_solver:Ch=%i:Iter=%i:Ant=%i"%(iChanSol,LMIter,iAnt),
                                       self._worker_solver,
                                       args=args)#,serial=True)
                    # End Antenna loop
                # End ChanLoop
                T3.timeit("APP.runJob")

                SmoothInterpLevel=(self.GD["Smoother"]["SmoothEnable"] and self.GD["Smoother"]["SmoothInterpLevel"])
                if SmoothInterpLevel and SmoothInterpLevel.lower()=="inner":
                    SmoothInterpMode=self.GD["Smoother"]["SmoothInterpMode"]
                    nd=self.SM.NDir
                    for iAnt in SortedWListAntSolve:
                        for iDir in range(nd):
                            APP.APP.runJob("workerSmooth:Ant=%i,iDir=%i"%(iAnt,iDir),
                                           self._workerSmooth,
                                           args=(iAnt,iDir))#,serial=True)
                    T.timeit("[%i] Smooth"%LMIter)
                    APP.APP.awaitJobResults("workerSmooth*")
                    #stop
                        
                if self.GD["Solvers"]["Sequential"]==0 and not isLastLMIter:
                    continue

                

                workers_res=APP.APP.awaitJobResults("worker_solver*")#, progress=1)
                T3.timeit("Await")
                
                
                T.timeit("EndChanAntennaloop")
                rmsFromDataList=[]
                DTs=np.zeros((self.VS.MS.na,),np.float32)
                    
                for iThisResult,ThisResult in enumerate(workers_res):
                    T4=ClassTimeIt.ClassTimeIt("  iter[%i][%i]"%(LMIter,iThisResult))
                    T4.disable()
                    iAnt,iChanSol,P,rmsFromData,InfoNoise,DT,Err,HasSolved = ThisResult
                    # print(iAnt,iChanSol,self.G1[iChanSol,iAnt,:].flat[0])
                    self.G[iChanSol,iAnt,:]=self.G1[iChanSol,iAnt,:]
                    self.HasSolvedArray[iChanSol,iAnt]=HasSolved
                    T4.timeit("Copy 0")
                    #G=self.G1[iChanSol,iAnt,:]
                    # if iAnt==0 and iChanSol==0:
                    #     print("AAA",iChanSol,iAnt,str(self.G[iChanSol,iAnt,:]))
                    #     stop
                    if rmsFromData!=None:
                        rmsFromDataList.append(rmsFromData)
                        
                    # T.timeit("[%i,%i] get"%(LMIter,iAnt))
                    #"TIMING DIFFERS BETWEEN SINGLE AND PARALLEL_NCPU=1"
                    #stop
                    #T.timeit("result_queue.get()")
                    #ThisG[iChanSol,iAnt][:]=G[:]
                    #self.G[iChanSol,iAnt][:]=G[:]
                    if type(P)!=type(None):
                        #P.fill(0.1)
                        ThisP[iChanSol,iAnt,:]=P[:]

                    T4.timeit("Copy 1")
                    DTs[iAnt]=DT
                    kapa=InfoNoise["kapa"]

                    if isLastLMIter:
                        self.DicoStd[iAnt,iChanSol].append(InfoNoise["std"])
                        self.DicoMax[iAnt,iChanSol].append(InfoNoise["max"])
                        self.DicoKeepKapa[iAnt,iChanSol].append(InfoNoise["kapa"])
                        self.SolsArray_Stats[self.iCurrentSol][iChanSol,iAnt][0]=InfoNoise["std"]
                        self.SolsArray_Stats[self.iCurrentSol][iChanSol,iAnt][1]=InfoNoise["max"]
                        self.SolsArray_Stats[self.iCurrentSol][iChanSol,iAnt][2]=InfoNoise["kapa"]
                        self.SolsArray_Stats[self.iCurrentSol][iChanSol,iAnt][3]=self.rms
                        
                    T4.timeit("Copy 2")
                    iResult+=1
                    if (kapa is not None)&(LMIter==0):
                        if kapa==-1.:
                            if len(self.DicoKapa[iAnt,iChanSol])>0:
                                kapa=self.DicoKapa[iAnt,iChanSol][-1]
                            else:
                                kapa=1.

                    T4.timeit("Copy 3")
                    if self.SolverType=="KAFCA":
                        self.DicoKapa[iAnt,iChanSol].append(kapa)
                        dt=.5
                        TraceResidList=self.DicoKapa[iAnt,iChanSol]
                        x=np.arange(len(TraceResidList))
                        expW=np.exp(-x/dt)[::-1]
                        expW/=np.sum(expW)
                        kapaW=np.sum(expW*np.array(TraceResidList))
                        ThisQ[iChanSol,iAnt][:]=(kapaW)*self.Q_Init[iChanSol,iAnt][:]
                    T4.timeit("Copy 4")


                    # if self.GD["Solutions"]["SmoothInterpMode"] is not None:
                    #     self.G[iChanSol,:]/=np.abs(self.G[iChanSol,:])
                    # self.G[iChanSol,:].fill(1000.)
                    
                    if self.SolverType=="KAFCA":
                        # self.P[iChanSol,iAnt,:]=ThisP[iChanSol,iAnt,:]
                        # self.Q[iChanSol,iAnt,:]=ThisQ[iChanSol,iAnt,:]
                        self.P[iChanSol,iAnt,:]=np.diag(np.diag(ThisP[iChanSol,iAnt,:]))
                        self.Q[iChanSol,iAnt,:]=np.diag(np.diag(ThisQ[iChanSol,iAnt,:]))
                        # nf,na,nd,nd=self.Q.shape
                        # P=np.zeros_like(self.P)
                        # Q=np.zeros_like(self.Q)
                        # # for iChan in range(nf):
                        # #     for iAnt in range(na):
                        # #         P[iChan,iAnt,:,:]=np.diag(np.diag(self.P[iChan,iAnt,:,:]))
                        # #         Q[iChan,iAnt,:,:]=np.diag(np.diag(self.Q[iChan,iAnt,:,:]))
                        # for iAnt in range(na):
                        #     P[iChanSol,iAnt,:,:]=np.diag(np.diag(self.P[iChanSol,iAnt,:,:]))
                        #     Q[iChanSol,iAnt,:,:]=np.diag(np.diag(self.Q[iChanSol,iAnt,:,:]))
                        # self.P[iChanSol,:]=P[iChanSol,:]
                        # self.Q[iChanSol,:]=Q[iChanSol,:]
                    T4.timeit("Copy 5")

                    
                    # if Err is not None:
                    #     import pylab
                    #     pylab.clf()
                    #     pylab.subplot(1,2,1)
                    #     pylab.scatter(self.SM.ClusterCat.SumI,Err)
                    #     pylab.title(iAnt)
                    #     pylab.subplot(1,2,2)
                    #     pylab.scatter(np.arange(Err.size),Err)
                    #     pylab.draw()
                    #     pylab.show(block=False)
                    #     pylab.pause(0.5)
                                
                    #print("kapa",kapaW)
                    # print self.Q[iChanSol,iAnt]
                    # print self.Q_Init[iChanSol,iAnt]
                    # print

                    # self.Q[iAnt][:]=(kapaW)**2*self.Q_Init[iAnt][:]*1e6
                    # QQ=NpShared.FromShared("%sSharedCovariance_Q"%self.IdSharedMem)[iAnt]
                    # print self.Q[iAnt]-QQ[iAnt]
                            
                    #self.Q[iAnt]=self.Q_Init[iAnt]
                    
                    #print iAnt,kapa,kapaW
                    #sig=np.sqrt(np.abs(np.array([np.diag(self.P[i]) for i in [iAnt]]))).flatten()
                    #print sig
                # End loop over queue results
                T3.timeit("ReadRes")

                T.timeit("[%i] OneIter"%LMIter)
                if len(rmsFromDataList)>0:
                    self.rmsFromData=np.min(rmsFromDataList)
                iResult=0

                # SmoothInterpLevel=self.GD["Solutions"]["SmoothInterpLevel"]
                # if SmoothInterpLevel is not None and SmoothInterpLevel.lower()=="inner":
                #     self.SmoothSols()
                #     T.timeit("[%i] Smooth"%LMIter)
                
                # pylab.clf()
                # pylab.subplot(2,1,1)
                # pylab.plot(DTs)
                # pylab.subplot(2,1,2)
                # pylab.plot(meanW)
                # pylab.draw()
                # pylab.show(False)
                # pylab.pause(0.1)
    
                if self.DoPlot==1:
                    import pylab
                    pylab.figure(1)
                    AntPlot=np.arange(self.VS.MS.na)#np.array(ListAntSolve)
                    pylab.clf()
                    pylab.plot(np.abs(ThisG[iChanSol,AntPlot].flatten()))
                    pylab.plot(np.abs(Gold[iChanSol,AntPlot].flatten()))
                        
                    if self.SolverType=="KAFCA":
                    
                        sig=[]
                        for iiAnt in AntPlot:
                            xx=np.array([np.diag(ThisP[iChanSol,iiAnt]) ])
                            if iiAnt in ListAntSolve:
                                sig.append(np.sqrt(np.abs(xx)).flatten().tolist())
                            else:
                                sig.append(np.zeros((xx.size,),ThisP.dtype).tolist())
                            
                        sig=np.array(sig).flatten()
                        pylab.plot(np.abs(ThisG[iChanSol,AntPlot].flatten())+sig,color="black",ls="--")
                        pylab.plot(np.abs(ThisG[iChanSol,AntPlot].flatten())-sig,color="black",ls="--")
                    pylab.title("Channel=%i"%iChanSol)
                    #pylab.ylim(0,2)
                    pylab.draw()
                    pylab.show(block=False)
                    pylab.pause(0.1)
                        
    
                T.timeit("[%i] Plot"%LMIter)

            # end Niter
            


            # #print self.P.ravel()
            # #if NDone==1: stop
            SmoothInterpLevel=(self.GD["Smoother"]["SmoothEnable"] and self.GD["Smoother"]["SmoothInterpLevel"])
            if SmoothInterpLevel and SmoothInterpLevel.lower()=="outer":
                SmoothInterpMode=self.GD["Smoother"]["SmoothInterpMode"]
                nd=self.SM.NDir
                for iAnt in SortedWListAntSolve:
                    for iDir in range(nd):
                        APP.APP.runJob("workerSmooth:Ant=%i,iDir=%i"%(iAnt,iDir),
                                       self._workerSmooth,
                                       args=(iAnt,iDir))#,serial=True)
                T.timeit("[%i] Smooth"%LMIter)
                APP.APP.awaitJobResults("workerSmooth*")
                # #stop
                # self.SmoothSols()
                T.timeit("Smooth")
                
            C0=(self.GD["VisData"]["FreePredictGainColName"]!=None or self.GD["SkyModel"]["FreeFullSub"])
            C1=(self.GD["VisData"]["FreePredictColName"]!=None)
            if C0 or C1:
                for iChanSol in range(self.VS.NChanJones):
                    for iAnt in SortedWListAntSolve:
                        DoFullPredict=True
                        
                        args=(iAnt,iChanSol,
                              DoCalcEvP[iChanSol],
                              tm,
                              self.rms,
                              DoEvP[iChanSol],
                              DoFullPredict,
                              (self.VS.iSliceSolveThis,self.VS.iChunkThis),
                              iiCount,LMIter)
                        APP.APP.runJob("worker_FreePredict:Ch=%i:Iter=%i:Ant=%i"%(iChanSol,LMIter,iAnt),
                                       self._worker_solver,
                                       args=args)#,serial=True)
                    
                APP.APP.awaitJobResults("worker_FreePredict*")
                    
            self.AppendGToSolArray()
            T.timeit("AppendGToSolArray")

            self.iCurrentSol+=1
            
            #_T=ClassTimeIt.ClassTimeIt("Plot")
            #_T.timeit()
            if self.DoPlot==2:
                S=self.GiveSols()
                #print S.G[-1,0,:,0,0]
                for ii in range(S.G.shape[1]):
                    self.Graph.subplot(ii)
                    self.Graph.imshow(np.abs(S.G[:,ii,:,0,0]).T)
                    #self.Graph.imshow(np.random.randn(*(S.G[:,ii,:,0,0]).shape))
                    self.Graph.text(0,0,self.VS.MS.StationNames[ii])
                self.Graph.draw()
                self.Graph.savefig()
            #_T.timeit()

            NDone,nt=self.pBarProgress
            intPercent=int(100*  NDone / float(nt))
            self.pBAR.render(NDone,nt)
            T.timeit("Ending")
            #if NDone==1:
            #    break
            self.VS.delCurrentSliceRelatedData()
            if self.VS_PredictCol is not None:
                self.VS_PredictCol.delCurrentSliceRelatedData()
            T.timeit("Delete")
            T2.timeit("One Slice Iter")


            
        # end while chunk

        # # To make the package more robust against memory leaks, we restart the worker processes every now and then.
        # # As a rule of thumb, we do this every major cycle, but no more often than every N minutes.
        # if time.time() > self.restart_time + 1200:
        #     APP.APP.restartWorkers()
        #     self.restart_time = time.time()
        
        APP.APP.restartWorkers()

        
        # if self.SolverType=="KAFCA":
        #     np.save("P.%s.npy"%self.GD["Solutions"]['OutSolsName'],np.array(self.PListKeep))
        #     np.savez("P.%s.npz"%self.GD["Solutions"]['OutSolsName'],
        #              ListP=np.array(self.PListKeep),
        #              ListQ=np.array(self.QListKeep))

            
        return True

    # ###############################################################################################################
    # ###############################################################################################################
    # ###############################################################################################################

    def setPredictParallel(self,PM):
        self.PM_Parallel=PM
        if self.GD["ImageSkyModel"]["BaseImageName"]!="" and self.GD["SkyModel"]["SkyModelCol"] is None and self.GD["FITSSkyModel"]["FITSSkyModel"] is None:
            self.PM_Parallel.InitGM(self.SM)

    def doCurrentFullPredict(self):
        DicoData=shared_dict.attach("DicoData_%i"%self.VS.iSliceSolveThis)
        
        ApplyTimeJones=None
        DicoPreApplyJones=shared_dict.attach("DicoPreApplyJones_Chunk%i"%self.VS.iChunkThis)
        if len(DicoPreApplyJones)>0:
            ApplyTimeJones=DicoPreApplyJones

        
        # DicoPredictData=shared_dict.attach("DicoPredictData")
        # DicoPredictData["PredictArray"]=np.zeros_like(DicoData["data"])
        
        PredictData=self.PM_Parallel.predictKernelPolCluster(DicoData,
                                                             self.SM,
                                                             iDirection=list(range(self.SM.NDir)),
                                                             ApplyTimeJones=ApplyTimeJones,
                                                             Progress=False)
        #DicoData["PredictDataDirs"]=PredictData
        #shared_dict.delDict("DicoPredictData")
    
    def _workerSmooth(self,iAnt,iDir):
        
        DicoCurrentGains=shared_dict.attach("DicoCurrentGains")
        DicoCurrentGains.reload()
        FreqsDone=DicoCurrentGains["FreqsDone"]

        if ("ATeam" in self.SM.ClusterCat.Name[iDir].decode("ascii")) and not self.GD["Smoother"]["AlsoSmoothATeam"]:
            return
        
        while not np.all(FreqsDone[iAnt,:]):
            # print("WAIT FOR iAnt=%i"%iAnt)
            time.sleep(0.1)
        
        G1=DicoCurrentGains["G1"]
        nChan,na,nd,nx,ny=G1.shape
        #if nx!=1 or ny!=1: stop
        for iPol in range(nx):
            for jPol in range(ny):
                x=G1[:,iAnt,iDir,iPol,jPol]
                nu=np.mean(self.VS.SolsFreqDomains,axis=1)
                SmoothInterpMode=self.GD["Smoother"]["SmoothInterpMode"]
                # x0=x.flatten().copy()
                
                xt=self.DicoSolsArray["G"][:,:,iAnt,iDir,iPol,jPol] # (NSols,nChan,na,nd,npolx,npoly)
                nDone=np.count_nonzero(self.DicoSolsArray["done"])
                n0=np.max([0,nDone-10])
                n1=nDone+1
                xts=xt[n0:n1].copy()
                xts[-1,:]=x
            
                r=ClassFitMultiModel.doFit_2D(xts.copy(),nu,Mode=SmoothInterpMode)#self.["Amp:Unit","Phase:TEC+CPhase"])
                #r=ClassFit_1D.doFit_1D(x.copy(),nu,Mode=SmoothInterpMode)#self.["Amp:Unit","Phase:TEC+CPhase"])
                G1[:,iAnt,iDir,iPol,jPol]=r[:]
        
                # x1=r.flatten().copy()
                
                # if iAnt==0:
                #     import pylab
                #     pylab.clf()
                #     pylab.subplot(1,2,1)
                #     pylab.plot(np.abs(x0),ls=":")
                #     pylab.plot(np.abs(x),ls="-")
                #     pylab.subplot(1,2,2)
                #     pylab.plot(np.angle(x0),ls=":")
                #     pylab.plot(np.angle(x),ls="-")
                #     pylab.draw()
                #     pylab.show(block=False)
                #     pylab.pause(0.1)



    def _worker_initJM(self,iAnt,
                       rms,
                       iSliceChunk):
        if self.SolverType=="CohJones":
            SolverClass=ClassSolverLM.ClassSolverLM
        elif self.SolverType=="KAFCA":
            SolverClass=ClassSolverEKF.ClassSolverEKF
        iChanSol=0
        LMIter=0
        iiCount=0

        iSlice=iSliceChunk[0]
        while not self.VS.DicoSemaphores["SemLoadSlice"][iSlice]:
            time.sleep(0.1)
        if self.VS_PredictCol is not None:
            while not self.VS_PredictCol.DicoSemaphores["SemLoadSlice"][iSlice]:
                time.sleep(0.1)

        # print("Init %i %i"%(iAnt,iSlice))
        ch0,ch1=self.VS.SolsToVisChanMapping[iChanSol]
        JM=SolverClass(self.SM,iAnt,PolMode=self.PolMode,
                       PM=self.PM,
                       PM_Compress=self.PM_Compress,
                       SM_Compress=self.SM_Compress,
                       IdSharedMem=self.IdSharedMem,
                       GD=self.GD,
                       ChanSel=(ch0,ch1),GlobalIter=(iiCount,LMIter),
                       iSliceChunk=iSliceChunk,
                       **dict(self.ConfigJacobianAntenna))
        JM.setDATA_Shared()
        JM.CalcKernelMatrix(rms)
        JM.AntennaDicoData_to_shm()
        self.PredictIsInit[iSlice,iAnt]=1
        # print("Init %i %i: Done"%(iAnt,iSlice))

        
    def _worker_solver(self,iAnt,iChanSol,DoCalcEvP,ThisTime,rms,DoEvP,DoFullPredict,iSliceChunk,iiCount,LMIter):
        
        if self.SolverType=="CohJones":
            SolverClass=ClassSolverLM.ClassSolverLM
        elif self.SolverType=="KAFCA":
            SolverClass=ClassSolverEKF.ClassSolverEKF
        T=ClassTimeIt.ClassTimeIt("     Worker iAnt,iChanSol,LMIter=%i,%i,%i"%(iAnt,iChanSol,LMIter))
        if not ((iAnt==0) and (iChanSol==0)):
            T.disable()
        T.disable()
        
        # if not (iAnt==0):
        #     T.disable()

        Mode="Solve"
        if DoFullPredict:
            Mode="DoFullPredict"

        
        ch0,ch1=self.VS.SolsToVisChanMapping[iChanSol]
        JM=SolverClass(self.SM,iAnt,PolMode=self.PolMode,
                       PM=self.PM,
                       PM_Compress=self.PM_Compress,
                       SM_Compress=self.SM_Compress,
                       IdSharedMem=self.IdSharedMem,
                       GD=self.GD,
                       ChanSel=(ch0,ch1),GlobalIter=(iiCount,LMIter),
                       iSliceChunk=iSliceChunk,
                       **dict(self.ConfigJacobianAntenna))

        
        iSlice=iSliceChunk[0]

        # #######################################
        JM.setDATA_Shared()
        if iChanSol==0 and LMIter==0 and not self.GD["Solvers"]["AsyncPredict"]:
            #print("INIT %i"%(iAnt))
            JM.CalcKernelMatrix(rms)
            self.PredictIsInit[iSlice,iAnt]=1
            #print("DONE %i"%(iAnt))
        else:
            while not self.PredictIsInit[iSlice,iAnt]:
                # print("WAIT %i %i %i"%(iAnt,iChanSol,LMIter))
                # print("WAIT %i %i %i"%(iAnt,iChanSol,LMIter))
                # print("WAIT %i %i %i"%(iAnt,iChanSol,LMIter))
                time.sleep(0.2)
            #print("---> OK %i"%(iAnt))
            JM.CalcKernelMatrix(rms)
        # # #######################################
        # # Initialise the predict for all channels/iChanSols
        # while not self.PredictIsInit[iSlice,iAnt]:
        #     time.sleep(0.2)
        # JM.setDATA_Shared()
        # JM.CalcKernelMatrix(rms)
        # # #################################
            
        JM.AntennaDicoData_to_shm()
        JM.SelectChannelKernelMat()

        T.timeit("setDATA_Shared")

        #shared_dict.attach("DicoSolsArray")
        DicoCurrentGains=shared_dict.attach("DicoCurrentGains")
        #DicoCurrentGains.reload()
        G=DicoCurrentGains["G"]
        

        GPrevious=G#NpShared.GiveArray("%sSharedGainsPrevious"%self.IdSharedMem)

        G0Iter=DicoCurrentGains["G0Iter"]

        #Q=NpShared.GiveArray("%sSharedCovariance_Q"%self.IdSharedMem)
        T.timeit("GiveArray")

        if self.SolverType=="CohJones":
            if Mode=="Solve":
                x,_,InfoNoise,HasSolved=JM.doLMStep(G[iChanSol],rms)
                T.timeit("LM")
            elif Mode=="DoFullPredict" and not JM.DoCompress: 
                # print("!!!!!!!!!!!!!!!!!!!")
                x=G[iChanSol,iAnt][:].copy()
                #Gc0=G.copy()
                # Gc0[iChanSol,iAnt][:]=x[:]
                Gc=G.copy()
                
                # Gc0.fill(2.)
                # NoZeroD=5
                # Gc.fill(0)
                # Gc[:,:,NoZeroD,:,:]=Gc0[:,:,NoZeroD,:,:]
                
                InfoNoise,HasSolved=None,None
                #Gc.fill(1)
                

                
                JM.PredictOrigFormat(Gc[iChanSol])
                T.timeit("FullPredict")
            else:
                stop
                
            Answer=[iAnt,
                    iChanSol,
                    None,
                    None,
                    InfoNoise,
                    0.,
                    None,
                    HasSolved]

        elif self.SolverType=="KAFCA":
            T0=time.time()
            DicoCov=shared_dict.attach("DicoCov")
            #DicoCov.reload()
            P=DicoCov["P"]
            evP=DicoCov["evP"]
            PPrevious=P#NpShared.GiveArray("%sSharedPPrevious"%self.IdSharedMem)
            #T.disable()
            # print("DoCalcEvP",DoCalcEvP)
            if DoCalcEvP:
                evP[iChanSol,iAnt]=JM.CalcMatrixEvolveCov(GPrevious[iChanSol],PPrevious[iChanSol],rms)
                
                #if evP[iChanSol,iAnt].max()>2: stop
                
                # if iAnt==51:
                #     M=(evP[iChanSol,iAnt]!=0)
                #     x=evP[iChanSol,iAnt][M].ravel()
                #     print("EVVV",x)
                #     print("EVVV",x)
                #     print("EVVV",x)
                T.timeit("Estimate Evolve")

            # EM=ClassModelEvolution(iAnt,
            #                        StepStart=3,
            #                        WeigthScale=2,
            #                        DoEvolve=False,
            #                        order=1,
            #                        sigQ=0.01)
            
            EM=ClassModelEvolution(iAnt,iChanSol,
                                   StepStart=0,
                                   WeigthScale=0.5,
                                   DoEvolve=True,
                                   BufferNPoints=10,
                                   sigQ=0.01,IdSharedMem=self.IdSharedMem)
            T.timeit("Init EM")

            Pa=None

            # Ga,Pa=EM.Evolve0(G,P,self.ThisTime)
            # if Ga!=None:
            #     G[iAnt]=Ga
            #     P[iAnt]=Pa
            
            # ThisP=P[iChanSol].copy()
            # ThisP.fill(0)
            # # print
            # # print ThisP.shape
            # # print
            # na,nd,_=ThisP.shape
            # for iAnt in range(na):
            #     for iDir in range(nd):
            #         ThisP[iAnt,iDir,iDir]=P[iChanSol][iAnt,iDir,iDir]
            # x,Pout,InfoNoise=JM.doEKFStep(G[iChanSol],ThisP,evP[iChanSol],rms,Gains0Iter=G0Iter)
            
            # if iAnt==51:
            #     print("KKKKK",iAnt,G[iChanSol].max(),P[iChanSol].max(),evP[iChanSol].max())
            
            if Mode=="Solve":
                x,Pout,InfoNoise,HasSolved=JM.doEKFStep(G[iChanSol],P[iChanSol],evP[iChanSol],rms,Gains0Iter=G0Iter)
                
                #print("[%i, %i] %f %f %f %f)"%(iAnt,iChanSol,np.max(x),np.max(G[iChanSol]),np.max(P[iChanSol]),np.max(evP[iChanSol])))
                #if HasSolved: stop
                #time.sleep(0.5)
                #if np.count_nonzero(np.abs(x)>10)>0: stop
                Err=None#JM.giveErrAntDir()
                rmsFromData=JM.rmsFromData
            
                T.timeit("EKFStep")
            elif Mode=="DoFullPredict" and not JM.DoCompress:
                x=G[iChanSol,iAnt][:].copy()
                JM.PredictOrigFormat(G[iChanSol])
                InfoNoise=None
                HasSolved=False
                Pout=None
                Err=None
                rmsFromData=None
                T.timeit("PredictOrigFormat")
            else:
                stop

            if DoEvP and HasSolved:
                #Pout0=Pout#.copy()
                Pa=EM.Evolve0(x,Pout)#,kapa=kapa)
                
                # if iAnt==51:
                #     print("###",iAnt,x.max(),Pa.max(),Pout0.max())
                
                T.timeit("Evolve")
            else:
                Pa=P[iChanSol,iAnt].copy()
                    
            #_,Pa=EM.Evolve(x,Pout,ThisTime)

            if type(Pa)!=type(None):
                Pout=Pa

            DT=time.time()-T0
            # L=[iAnt,
            #    iChanSol,
            #    x,
            #    Pout,
            #    rmsFromData,
            #    InfoNoise,
            #    DT]
            # self.result_queue.put(L)
                
            Answer= [iAnt,
                     iChanSol,
                     Pout,
                     rmsFromData,
                     InfoNoise,
                     DT,
                     Err,
                     HasSolved]

            
        if self.GD["Solvers"]["Sequential"]==1:
            self.DicoCurrentGains["G1"][iChanSol,iAnt].flat[:]=x.flat[:]
        else:
            self.DicoCurrentGains["G"][iChanSol,iAnt].flat[:]=x.flat[:]

        # if self.GD["Solvers"]["SmoothInterpLevel"] is None:
        #     self.DicoCurrentGains["G"][iChanSol,iAnt].flat[:]=x.flat[:]
        #     self.DicoCurrentGains["G1"][iChanSol,iAnt].flat[:]=x.flat[:]

        #print("done: ",iAnt,iChanSol)
        self.DicoCurrentGains["FreqsDone"][iAnt,iChanSol]=1
        # if iChanSol==0 and iAnt==0:
        #     print("BBB",iChanSol,iAnt,x)
        #     print("BBB1",iChanSol,iAnt,self.DicoCurrentGains["G1"][iChanSol,iAnt].flat[:])
        return Answer
    


        
