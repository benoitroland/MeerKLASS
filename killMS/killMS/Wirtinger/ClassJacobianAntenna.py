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
# from killMS.Array import NpShared
from killMS.Predict.PredictGaussPoints_NumExpr5 import ClassPredict
import os
from killMS.Data import ClassVisServer
#from Sky import ClassSM
from killMS.Array import ModLinAlg
from killMS.Other import ClassTimeIt
#from killMS.Array.Dot import NpDotSSE
from . import ClassAverageMachine
from SkyModel.Sky import ClassSM
from DDFacet.Other import logger
log=logger.getLogger("ClassJacobianAntenna")
from DDFacet.Array import shared_dict
#from DDFacet.Other.AsyncProcessPool import APP

def testLM():
    import pylab
    SM=ClassSM.ClassSM("../TEST/ModelRandom00.txt.npy")
    rabeam,decbeam=SM.ClusterCat.ra,SM.ClusterCat.dec
    LofarBeam=("AE",5,rabeam,decbeam)
    VS=ClassVisServer.ClassVisServer("../TEST/0000.MS/")#,LofarBeam=LofarBeam)
    MS=VS.MS
    SM.Calc_LM(MS.rac,MS.decc)



    nd=SM.NDir
    npol=4
    na=MS.na
    Gains=np.zeros((na,nd,npol),dtype=np.complex64)
    Gains[:,:,0]=1
    Gains[:,:,-1]=1
    #Gains+=(np.random.randn(*Gains.shape)*0.5+1j*np.random.randn(*Gains.shape)
    #Gains=np.random.randn(*Gains.shape)+1j*np.random.randn(*Gains.shape)
    #GainsOrig=Gains.copy()
    ###############
    #GainsOrig=np.load("Rand.npz")["GainsOrig"]
    #Gains*=1e-3
    #Gains[:,0,:]=GainsOrig[:,0,:]
    ###############

    PolMode="IFull"
    # ### Scalar gains
    # PolMode="Scalar"
    # g=np.random.randn(*(Gains[:,:,0].shape))+1j*np.random.randn(*(Gains[:,:,0].shape))
    # g=g.reshape((na,nd,1))
    # Gains*=g
    # ####
    
    #Gains[:,2,:]=0
    #Gains[:,1,:]=0
    
    #Gains[:,1::,:]=0
    


    DATA=VS.GiveNextVis()


    # Apply Jones
    PM=ClassPredict(Precision="S")
    #DATA["data"]=PM.predictKernelPolCluster(DATA,SM,ApplyJones=Gains)
    
    ############################
    iAnt=0
    JM=ClassJacobianAntenna(SM,iAnt,PolMode=PolMode)
    JM.setDATA(DATA)
    JM.CalcKernelMatrix()

    if PolMode=="Scalar":
        Gains=Gains[:,:,0].reshape((na,nd,1))
        G=Gains.copy().reshape((na,nd,1,1))
    else:
        G=Gains.copy().reshape((na,nd,2,2))
        

    y=JM.GiveDataVec()
    xtrue=JM.GiveSubVecGainAnt(Gains).flatten()
    x=Gains
    #################
    Radd=np.random.randn(*(G[iAnt].shape))#*0.3
    #np.savez("Rand",Radd=Radd,GainsOrig=GainsOrig)
    #################
    #Radd=np.load("Rand.npz")["Radd"]

    G[iAnt]+=Radd


    print("start")
    for i in range(10):
        xbef=G[iAnt].copy()
        x=JM.doLMStep(G)
        G[iAnt]=x
        
        pylab.figure(1)
        pylab.clf()
        pylab.plot(np.abs(xtrue.flatten()))
        pylab.plot(np.abs(x.flatten()))
        pylab.plot(np.abs(xtrue.flatten())-np.abs(x.flatten()))
        pylab.plot(np.abs(xbef.flatten()))
        pylab.ylim(-2,2)
        pylab.draw()
        pylab.show(False)

        # pylab.figure(3)
        # pylab.clf()
        # pylab.subplot(1,3,1)
        # pylab.imshow(np.abs(JM.Jacob)[0:20],interpolation="nearest")
        # pylab.subplot(1,3,2)
        # pylab.imshow(np.abs(JM.JHJ),interpolation="nearest")
        # pylab.subplot(1,3,3)
        # pylab.imshow(np.abs(JM.JHJinv),interpolation="nearest")
        # pylab.draw()
        # pylab.show(False)

    stop    
    



class ClassJacobianAntenna():
    def __init__(self,SM,iAnt,PolMode="IFull",Precision="S",PrecisionDot="D",IdSharedMem="",
                 PM=None,
                 PM_Compress=None,
                 SM_Compress=None,
                 GD=None,NChanSols=1,ChanSel=None,
                 GlobalIter=(0,0),
                 iSliceChunk=None,
                 VisToSolsChanMapping=None,
                 **kwargs):
        T=ClassTimeIt.ClassTimeIt("  InitClassJacobianAntenna")
        T.disable()
        self.VisToSolsChanMapping=VisToSolsChanMapping
        self.NChanSol=np.unique(self.VisToSolsChanMapping).size
        self.maxFlagFrac=0.2
        self.Err=None
        self.iiCount,self.LMIter=GlobalIter
        self.iCurrentVisTime,self.iTimeChunk=iSliceChunk
        self.IndicesSel0=None
        self.IndicesSel1=None
        self.IndicesSel2=None
        self.IndicesSel3=None
        self.DATA_ColumnPredict=None        
        self.ChanSel=ChanSel
        self.GD=GD
        self.IdSharedMem=IdSharedMem
        self.PolMode=PolMode
        #self.PM=ClassPredict(Precision="S")
        self.Rinv_flat=None
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
        
        self.PM=PM
        self.SM=SM
        T.timeit("Init0")
        if PM==None:
            self.PM=ClassPredict(Precision=Precision,
                                 DoSmearing=self.DoSmearing,
                                 IdMemShared=IdSharedMem,
                                 BeamAtFacet=(self.GD["Beam"]["BeamAt"].lower() == "facet"))
            
            if self.GD["ImageSkyModel"]["BaseImageName"]!="":
                self.PM.InitGM(self.SM)

        self.PM_Compress=PM_Compress
        self.SM_Compress=SM_Compress
        self.AverageMachine=None
        self.DoCompress=False
        self.DoMergeStations = ((self.GD["Compression"]["MergeStations"] is not None) and (self.GD["Compression"]["MergeStations"]!=""))
        if self.SM_Compress or self.DoMergeStations:
            self.AverageMachine=ClassAverageMachine.ClassAverageMachine(self.GD,
                                                                        self.PM_Compress,
                                                                        self.SM_Compress,
                                                                        DicoMergeStations=self.DicoMergeStations,
                                                                        VisToSolsChanMapping=VisToSolsChanMapping)
            self.DoCompress=True
            self.NDirAvg = self.SM_Compress.NDir
        
        T.timeit("PM")
        if PrecisionDot=="D":
            self.CType=np.complex128
            self.FType=np.float64

        if PrecisionDot=="S":
            self.CType=np.complex64
            self.FType=np.float32

        self.CType=np.complex128
        self.TypeDot="Numpy"
        # self.TypeDot="SSE"

        self.iAnt=int(iAnt)
        self.KernelSharedName="KernelMat.A_%i.Slice_%i"%(self.iAnt,self.iCurrentVisTime)
        self.DicoMainDataName="DicoData_%i"%self.iCurrentVisTime
        self.DicoColPredictDataName="ColumnPredictVS_DicoData_%i"%self.iCurrentVisTime
        
        self.NameDicoPreApplyJones="DicoPreApplyJones_Chunk%i"%self.iTimeChunk
        self.DicoPredictedVisName="DicoPredictedVis_Chunk%i"%self.iTimeChunk
        self.DicoThisDataChunkName="DataChunk_Chunk%i"%self.iTimeChunk
        
        self.NChanSols=NChanSols
        
        if self.PolMode=="IFull":
            self.NJacobBlocks_X=2
            self.NJacobBlocks_Y=2
            self.npolData=4
            
        elif self.PolMode=="Scalar":
            self.NJacobBlocks_X=1
            self.NJacobBlocks_Y=1
            self.npolData=1
        
        elif self.PolMode=="IDiag":
            self.NJacobBlocks_X=2
            self.NJacobBlocks_Y=1
            self.npolData=2
        

        self.Reinit()
        T.timeit("rest")

    def Reinit(self):
        self.HasKernelMatrix=False
        self.LQxInv=None

    def GiveSubVecGainAnt(self,GainsIn):
        # if (GainsIn.size==self.NDir*2*2): return GainsIn.copy()
        Gains=GainsIn.copy().reshape((self.na,self.NDir,self.NJacobBlocks_X,self.NJacobBlocks_Y))[self.iAnt]
        return Gains
        
    def setDATA(self,DATA):
        self.DATA=DATA
        
    def setDATA_Shared(self):
        # SharedNames=["SharedVis.freqs","SharedVis.times","SharedVis.A1","SharedVis.A0","SharedVis.flags","SharedVis.infos","SharedVis.uvw","SharedVis.data"]
        # self.DATA={}
        # for SharedName in SharedNames:
        #     key=SharedNames.split(".")[1]
        #     self.DATA[key]=NpShared.GiveArray(SharedName)

        T=ClassTimeIt.ClassTimeIt("  setDATA_Shared")
        T.disable()
        #self.DATA=NpShared.SharedToDico("%sSharedVis"%self.IdSharedMem)
        #print("JCOB ATTACH 0")
        self.DATA=shared_dict.attach(self.DicoMainDataName)

        if self.SM.Type=="Column":
            #print("JCOB ATTACH 1 ")
            #os.system("ls /dev/shm/ddf.*/CulumnPredictVS_DicoData_0")
            self.DATA_ColumnPredict=shared_dict.attach(self.DicoColPredictDataName,readwrite=False)



        _,self.NChanMS,_=self.DATA["data"].shape
        if self.ChanSel==None:
            self.ch0=0
            self.ch1=self.NChanMS
        else:
            self.ch0,self.ch1=self.ChanSel

        self.SharedDataDicoName="DicoData.A_%i.ch_%i-%i.Slice_%i"%(self.iAnt,self.ch0,self.ch1,self.iCurrentVisTime)
        
        self.NChanData=self.ch1-self.ch0

        T.timeit("SharedToDico0")
        #DicoBeam=NpShared.SharedToDico("%sPreApplyJones"%self.IdSharedMem)

        # DicoPreApplyJones=shared_dict.attach("DicoPreApplyJones")
        # #DicoBeam=DicoThisDataChunk.get("PreApplyJones",None)
        # #DicoBeam=shared_dict.attach(self.NameDicoPreApplyJones)
        # T.timeit("SharedToDico1")
        # if len(DicoPreApplyJones)>0:
        #     self.DATA["DicoPreApplyJones"]=DicoBeam
        #     self.DATA["DicoClusterDirs"]=shared_dict.attach("DicoClusterDirs")


        T.timeit("SharedToDico2")


        #self.DATA["UVW_RefAnt"]=NpShared.GiveArray("%sUVW_RefAnt"%self.IdSharedMem)



           
            

                                        
    def JHJinv_x(self,Gains):
        G=[]
        #nd,_,_=Gains.shape
        Gains=Gains.reshape((self.NDir,self.NJacobBlocks_X,self.NJacobBlocks_Y))
        for polIndex in range(self.NJacobBlocks_X):
            Gain=Gains[:,polIndex,:]
            #print("JHJinv_x: %i %s . %s "%(polIndex,str(self.L_JHJinv[polIndex].shape),str(Gain.flatten().shape)))
            Vec=np.dot(self.L_JHJinv[polIndex],Gain.flatten())
            Vec=Vec.reshape((self.NDir,1,self.NJacobBlocks_Y))
            G.append(Vec)
            
        Gout=np.concatenate(G,axis=1)
        #print("JHJinv_x: Gout %s "%(str(Gout.shape)))
        
        return Gout.flatten()



    def Msq_x(self,LM,Gains):
        G=[]
        Gains=Gains.reshape((self.NDir,self.NJacobBlocks_X,self.NJacobBlocks_Y))
        for polIndex in range(self.NJacobBlocks_X):
            Gain=Gains[:,polIndex,:]
            #print("Msq_x: %i %s . %s"%(polIndex,str(LM[polIndex].shape),str(Gain.flatten().shape)))
            Vec=np.dot(LM[polIndex],Gain.flatten())
            Vec=Vec.reshape((self.NDir,1,self.NJacobBlocks_Y))
            G.append(Vec)
            
        Gout=np.concatenate(G,axis=1)
        #print("Msq_x: Gout %s "%(str(Gout.shape)))
        
        return Gout.flatten()




    def JH_z(self,zin):
        #z=zin.reshape((self.NJacobBlocks,zin.size/self.NJacobBlocks))
        #z=zin.reshape((1,zin.size))
        Gains=np.zeros((self.NDir,self.NJacobBlocks_X,self.NJacobBlocks_Y),self.CType)

        if self.DoCompress:
            flags_key="flags_flat_avg"
            if self.DoMergeStations:
                flags_key="flags_flat_avg_merged"
        else:
            flags_key="flags_flat"
            

        for polIndex in range(self.NJacobBlocks_X):

            Jacob=self.LJacob[polIndex]

            
            flags=self.DicoData[flags_key][polIndex]
            ThisZ=zin[polIndex][flags==0]#self.DicoData["flags_flat"[polIndex]
            
            J=Jacob[flags==0]


            nrow,_=J.shape
            ThisZ=ThisZ.flatten()
            C_CohJones=(self.GD["Solvers"]["SolverType"]=="CohJones")
            if (self.Rinv_flat is not None) and C_CohJones:
                Rinv=self.Rinv_flat[polIndex][flags==0].reshape((nrow,1))
                if self.TypeDot=="Numpy":
                    Gain=np.dot(J.T.conj(),Rinv.flatten()*ThisZ)
                elif self.TypeDot=="SSE":
                    ThisZ=ThisZ.reshape((1,ThisZ.size))
                    ThisZ.flat[:]*=Rinv.flat[:]
                    JTc=self.LJacobTc[polIndex]#.copy()
                    Gain=NpDotSSE.dot_A_BT(JTc,ThisZ)
            else:
                if self.TypeDot=="Numpy":
                    Gain=np.dot(J.T.conj(),ThisZ)
                elif self.TypeDot=="SSE":
                    ThisZ=ThisZ.reshape((1,ThisZ.size))
                    JTc=self.LJacobTc[polIndex]#.copy()
                    Gain=NpDotSSE.dot_A_BT(JTc,ThisZ)


            Gains[:,polIndex,:]=Gain.reshape((self.NDir,self.NJacobBlocks_Y))

        return Gains

    # def GiveDataVec(self):
    #     y=[]
    #     yf=[]
    #     for polIndex in range(self.NJacobBlocks):
    #         DataVec=self.DicoData["data"][:,:,polIndex,:].flatten()
    #         Flags=self.DicoData["flags"][:,:,polIndex,:].flatten()

    #         y.append(DataVec)
    #         yf.append(Flags)

    #     y=np.concatenate(y)
    #     yf=np.concatenate(yf)
    #     return y,yf


    def J_x(self,Gains):
        z=[]
        Gains=Gains.reshape((self.NDir,self.NJacobBlocks_X,self.NJacobBlocks_Y))
        for polIndex in range(self.NJacobBlocks_X):
            Jacob=self.LJacob[polIndex]
            
            Gain=Gains[:,polIndex,:].flatten()

            #flags=self.DicoData["flags_flat"][polIndex]
            J=Jacob#[flags==0]
            # print(J.shape, Gain.shape)

            # # Numpy

            if self.TypeDot=="Numpy":
                Z=np.dot(J,Gain)
            elif self.TypeDot=="SSE":
                Gain=Gain.reshape((1,Gain.size))
                Z=NpDotSSE.dot_A_BT(J,Gain).ravel()

            z.append(Z)



        z=np.array(z)
        return z

    def PredictOrigFormat(self,GainsIn):
        if self.GD["VisData"]["FreePredictGainColName"]!=None or self.GD["SkyModel"]["FreeFullSub"]:
            self.PredictOrigFormat_Type(GainsIn,Type="Gains")
        if self.GD["VisData"]["FreePredictColName"]!=None:
            self.PredictOrigFormat_Type(GainsIn,Type="NoGains")


    def PredictOrigFormat_Type(self,GainsIn,Type="Gains"):
        T=ClassTimeIt.ClassTimeIt("     -- PredictOrigFormat ")
        T.disable()
        #print("    COMPUTE PredictOrigFormat")
        Gains=GainsIn.copy()
        na,nd,_,_=Gains.shape
        #Type="NoGains"

        DicoPredictedVis=shared_dict.attach(self.DicoPredictedVisName)
        
        OpCol=(lambda x:x)
        D0=np.zeros_like(self.DicoData["data_allpol"])
        if Type=="NoGains":
            if self.PolMode=="Scalar":
                Gains=np.ones((na,nd,1,1),np.complex64)
            elif self.PolMode=="IDiag":
                Gains=np.ones((na,nd,2,1),np.complex64)
            else:
                Gains=np.zeros((na,nd,2,2),np.complex64)
                Gains[:,:,0,0]=1
                Gains[:,:,1,1]=1
            NameShmData="PredictedData"
            NameShmIndices="IndicesData"
        elif Type=="Gains":
            DataName=None
            if ":" in self.GD["VisData"]["FreePredictGainColName"]: # for example --FreePredictGainColName=KMS_SUB_SCALAR:data-ATeam
                ColName,Op=self.GD["VisData"]["FreePredictGainColName"].split(":")
                if "-" in Op:
                    OpCol=(lambda x:-x)
                    DataName,DirName=Op.split("-")
                    if DataName!="data":
                        stop
                    D0=self.DicoData["data_allpol"]
                    DirsName=[self.SM.ClusterCat.Name[i].decode("ascii") for i in range(self.SM.ClusterCat.Name.size)]
                    for iDir in range(nd):
                        if DirName not in DirsName[iDir]:
                            Gains[:,iDir]=0
                else:
                    stop


            NameShmData="PredictedDataGains"
            NameShmIndices="IndicesDataGains"

        PredictedData=DicoPredictedVis[NameShmData]
        Indices=DicoPredictedVis[NameShmIndices]
        
        Ga=self.GiveSubVecGainAnt(Gains).copy()
        
        T.timeit("Init")
        self.CalcJacobianAntenna(Gains)
        T.timeit("CalcJacobianAntenna")
        #self.PrepareJHJ_LM()
        zp=self.J_x(Ga)#self.DicoData["data_flat"]#
        T.timeit("zp")
        
        z=self.DicoData["data_flat"]#self.GiveDataVec()
        f=(self.DicoData["flags_flat"]==0)
        zr=(z-zp)
        zrf=zr[f]
        if zrf.size>0:
            self.rmsFromData=np.std(zrf)
            
        DicoData=self.DicoData

        nr,nch,_,_=DicoData["flags"].shape
            
        indRowsThisChunk=self.DATA["indRowsThisChunk"]
        indOrig=DicoData["indOrig"]
        indThis=DicoData["indThis"]
        
        T.timeit("IndThis")

        if self.IndicesSel0 is None:
            # self.IndicesSel0=Indices[indRowsThisChunk,:,:][indOrig,self.ch0:self.ch1,0].ravel()
            self.IndicesSel0=Indices[indRowsThisChunk[indOrig],self.ch0:self.ch1,0].ravel()
            # stop
            T.timeit("IndicesSel0")
            self.IndicesSel3=self.IndicesSel0+3#Indices[indRowsThisChunk,:,:][indOrig,self.ch0:self.ch1,3].ravel()
            T.timeit("IndicesSel3")
            if self.PolMode=="IFull":
                self.IndicesSel1=self.IndicesSel0+1#Indices[indRowsThisChunk,:,:][indOrig,self.ch0:self.ch1,1].ravel()
                self.IndicesSel2=self.IndicesSel0+2#Indices[indRowsThisChunk,:,:][indOrig,self.ch0:self.ch1,2].ravel()
                
        T.timeit("IndicesSel")
        
        D=np.rollaxis(zp.reshape(self.NJacobBlocks_X,nr,nch,self.NJacobBlocks_Y),0,3).reshape(nr,nch,self.NJacobBlocks_X,self.NJacobBlocks_Y)
        T.timeit("roll")

        if self.PolMode=="Scalar":
            # PredictedData.ravel()[IndicesSel0]=D[indThis,:,0,0].ravel()
            PredictedData.flat[self.IndicesSel0]=D0[indThis,...,0,0].ravel()+OpCol(D[indThis,:,0,0].ravel())
            PredictedData.flat[self.IndicesSel3]=D0[indThis,...,1,1].ravel()+OpCol(D[indThis,:,0,0].ravel())
        elif self.PolMode=="IDiag":
            PredictedData.flat[self.IndicesSel0]=D0[indThis,...,0,0].ravel()+OpCol(D[indThis,:,0,0].ravel())
            PredictedData.flat[self.IndicesSel3]=D0[indThis,...,1,1].ravel()+OpCol(D[indThis,:,1,0].ravel())
        elif self.PolMode=="IFull":
            PredictedData.flat[self.IndicesSel0]=D0[indThis,...,0,0].ravel()+OpCol(D[indThis,:,0,0].ravel())
            PredictedData.flat[self.IndicesSel1]=D0[indThis,...,0,1].ravel()+OpCol(D[indThis,:,0,1].ravel())
            PredictedData.flat[self.IndicesSel2]=D0[indThis,...,1,0].ravel()+OpCol(D[indThis,:,1,0].ravel())
            PredictedData.flat[self.IndicesSel3]=D0[indThis,...,1,1].ravel()+OpCol(D[indThis,:,1,1].ravel())
        T.timeit("flat")

        
        
        # if Type=="NoGains":
        #     stop
        
        # d0=self.DATA["data"]#[indOrig,:,0]
        # #d1=D[indThis,:,0,0]
        # d2=PredictedData[indRowsThisChunk,:,:]#[indOrig,:,0]

        # pylab.clf()
        # pylab.plot(d0[:,2,0].real)
        # # pylab.plot(d1[:,2,0].real)
        # pylab.plot(d2[:,2,0].real)
        # pylab.plot((d0-d2)[:,2,0].real)
        # pylab.draw()
        # pylab.show(False)
        # pylab.pause(0.1)
        # # stop




    def CalcJacobianAntenna(self,GainsIn):
        if not(self.HasKernelMatrix): stop
        #if np.count_nonzero(np.isnan(GainsIn))>0: stop
        iAnt=self.iAnt
        NDir=self.NDir

        if self.DoCompress:
            n4vis=self.n4vis_Avg
        else:
            n4vis=self.n4vis
            
        #print("n4vis",n4vis)
        na=self.na
        #print(GainsIn.shape,na,NDir,self.NJacobBlocks,self.NJacobBlocks)
        Gains=GainsIn.reshape((na,NDir,self.NJacobBlocks_X,self.NJacobBlocks_Y))
        
        Jacob=np.zeros((n4vis,self.NJacobBlocks_Y,NDir,self.NJacobBlocks_Y),self.CType)
        Jacob.fill(0)
        
        if (self.PolMode=="IFull")|(self.PolMode=="Scalar"):
            self.LJacob=[Jacob]*self.NJacobBlocks_X
        elif self.PolMode=="IDiag":
            self.LJacob=[Jacob,Jacob.copy()]
        LJacob=self.LJacob

        if self.DoCompress:
            A1=self.DicoData["A1_Avg"]
        else:
            A1=self.DicoData["A1"]
        
        for iDir in range(NDir):
            G=Gains[A1,iDir].conj()
            
            K_XX=self.K_XX[iDir]
            K_YY=self.K_YY[iDir]

            nr=G.shape[0]

            if self.PolMode=="Scalar":
                J0=Jacob[:,0,iDir,0]
                g0_conj=G[:,0,0].reshape((nr,1))
                J0[:]=(g0_conj*K_XX).reshape((K_XX.size,))
            elif self.PolMode=="IFull":
                J0=Jacob[:,0,iDir,0]
                g0_conj=G[:,0,0].reshape((nr,1))
                J0[:]=(g0_conj*K_XX).reshape((K_XX.size,))
                J1=Jacob[:,0,iDir,1]
                J2=Jacob[:,1,iDir,0]
                J3=Jacob[:,1,iDir,1]
                g1_conj=G[:,1,0].reshape((nr,1))
                g2_conj=G[:,0,1].reshape((nr,1))
                g3_conj=G[:,1,1].reshape((nr,1))
                J1[:]=(g2_conj*K_YY).reshape((K_XX.size,))
                J2[:]=(g1_conj*K_XX).reshape((K_XX.size,))
                J3[:]=(g3_conj*K_YY).reshape((K_XX.size,))
            elif self.PolMode=="IDiag":
                J0=LJacob[0][:,0,iDir,0]
                g0_conj=G[:,0,0].reshape((nr,1))
                J0[:]=(g0_conj*K_XX).reshape((K_XX.size,))
                J1=LJacob[1][:,0,iDir,0]
                g1_conj=G[:,1,0].reshape((nr,1))
                J1[:]=(g1_conj*K_YY).reshape((K_XX.size,))


        for J in LJacob:
            J.shape=(n4vis*self.NJacobBlocks_Y,NDir*self.NJacobBlocks_Y)
            
        if self.DoCompress:
            flags_key="flags_flat_avg"
        else:
            flags_key="flags_flat"

        self.LJacobTc=[]
        for polIndex in range(self.NJacobBlocks_X):
            flags=self.DicoData[flags_key][polIndex]
            J=self.LJacob[polIndex][flags==0]
            self.LJacobTc.append(J.T.conj().copy())

        self.L_JHJ=[]
        for polIndex in range(self.NJacobBlocks_X):
            flags=self.DicoData[flags_key][polIndex]
            J=self.LJacob[polIndex][flags==0]
            nrow,_=J.shape
            self.nrow_nonflagged=nrow
            # JH=J.T.conj()
            JH=self.LJacobTc[polIndex]
            if self.Rinv_flat is not None:
                
                Rinv=self.Rinv_flat[polIndex][flags==0].reshape((nrow,1))
                
                if self.TypeDot=="Numpy":
                    JHJ=np.dot(JH,Rinv*J)
                elif self.TypeDot=="SSE":
                    RinvJ_T=(Rinv*J).T.copy()
                    JTc=self.LJacobTc[polIndex]#.copy()
                    JHJ=NpDotSSE.dot_A_BT(JTc,RinvJ_T)

                # if np.count_nonzero(flags==0)>0:
                #     L_JHJr=[]
                #     L_y=[]
                #     L_x=[]
                #     NTry=10
                #     l=self.SM.ClusterCat.l
                #     m=self.SM.ClusterCat.m
                #     u,v,w=self.DicoData["uvw"].T
                #     f=self.DicoData["flags"]
                #     freqs=self.DicoData["freqs"]
                    
                #     ClusterCat=self.SM.ClusterCat.view(np.recarray)
                #     SourceCat=np.zeros((ClusterCat.shape[0],),dtype=[('Name', 'S200'), ('ra', '<f8'), ('dec', '<f8'), ('Sref', '<f8'),
                #                                                      ('I', '<f8'), ('Q', '<f8'), ('U', '<f8'), ('V', '<f8'), ('RefFreq', '<f8'),
                #                                                      ('alpha', '<f8'), ('ESref', '<f8'), ('Ealpha', '<f8'), ('kill', '<i8'),
                #                                                      ('Cluster', '<i8'), ('Type', '<i8'), ('Gmin', '<f8'), ('Gmaj', '<f8'),
                #                                                      ('Gangle', '<f8'), ('Select', '<i8'), ('l', '<f8'), ('m', '<f8'), ('Exclude', '<i8')])
                #     SourceCat=SourceCat.view(np.recarray)
                #     SourceCat.ra[:]=ClusterCat.ra[:]
                #     SourceCat.dec[:]=ClusterCat.dec[:]
                #     SourceCat.RefFreq[:]=100.e6
                #     SourceCat.I[:]=1.
                #     SourceCat.Cluster=np.arange(ClusterCat.shape[0])
                #     np.save("SM_Compress.%i.npy"%self.iAnt,SourceCat)
                #     SM_Cluster=ClassSM.ClassSM("SM_Compress.%i.npy"%self.iAnt,DoPrint=False)
                #     SM_Cluster.Calc_LM(self.SM.rac,self.SM.decc)

                #     PM=ClassPredict(Precision="S",
                #                     IdMemShared=self.IdSharedMem)
                    
                #     L_k=[]
                #     for iDir in range(self.NDir):
                #         K=PM.predictKernelPolCluster(self.DicoData,SM_Cluster,iDirection=iDir)
                #         L_k.append(K[...,0].reshape((-1,))[flags==0])

                #     K=np.array(L_k).T
                #     Th=np.max(np.array([ClusterCat.SumI,np.ones((self.NDir,))*5]),axis=0)
                    
                #     for iTry in range(NTry):
                #         rd=(np.random.randn(self.NDir)+1j*np.random.randn(self.NDir))/np.sqrt(2)
                #         rd=rd.reshape((1,self.NDir))
                #         Jc=J.copy()*rd
                #         #Jc+=(K*rd)*Th.reshape((1,self.NDir))*1e-3#*self.DicoData["rms"]#1e-3
                #         JHJr=np.dot(Jc.T.conj(),Jc)
                #         JrHJ=np.dot(Jc.T.conj(),J)
                #         #JHJr.append(np.dot(Jc.T.conj(),Rinv*Jc))
                #         A=JHJr
                #         Ai=ModLinAlg.invSVD(A)
                #         g=Gains[self.iAnt,:]
                #         x=g.reshape((-1,1))+(np.random.randn(A.shape[0],1)+1j*np.random.randn(A.shape[0],1))/np.sqrt(2)*0.1
                #         y=Ai.dot(JrHJ).dot(x)
                #         L_x.append(x)
                #         L_y.append(y)
                #     self.Err=np.var(np.array(L_y)-np.array(L_x),axis=0)
             
            else:
                if self.TypeDot=="Numpy":
                    JHJ=np.dot(JH,J)
                elif self.TypeDot=="SSE":
                    J_T=J.T.copy()
                    JTc=self.LJacobTc[polIndex]#.copy()
                    JHJ=NpDotSSE.dot_A_BT(JTc,J_T)
                    
            if np.count_nonzero(np.isnan(JHJ))>0: stop
            
            self.L_JHJ.append(self.CType(JHJ))
            
        if self.DoMergeStations:
            self.LJacob=self.AverageMachine.MergeAntennaJacobian(self.DicoData,LJacob)
            #print(self.iAnt,len(self.LJacob))

        # self.JHJinv=np.linalg.inv(self.JHJ)
        # self.JHJinv=np.diag(np.diag(self.JHJinv))
        
    def CalcKernelMatrix(self,rms=0.):
        # Out[28]: ['freqs', 'times', 'A1', 'A0', 'flags', 'uvw', 'data']
        T=ClassTimeIt.ClassTimeIt("CalcKernelMatrix Ant=%i"%self.iAnt)
        T.disable()
        DATA=self.DATA
        iAnt=self.iAnt
        na=int(DATA['infos'][0])
        self.na=na
        NDir=self.SM.NDir
        self.NDir=NDir
        self.iAnt=iAnt

        T.timeit("stuff")
        
        self.DicoData=self.GiveData(DATA,iAnt,rms=rms)
       
        T.timeit("data")
        # self.Data=self.DicoData["data"]
        # self.A1=self.DicoData["A1"]
        # print("AntMax1",self.SharedDataDicoName,np.max(self.A1))
        # print(self.DicoData["A1"])
        # print("AntMax0",self.SharedDataDicoName,np.max(self.DicoData["A0"]))
        # print(self.DicoData["A0"])
        nrows,nchan,_,_=self.DicoData["flags"].shape
        n4vis=nrows*nchan
        self.n4vis=n4vis

        # if self.DoCompress:
        #     nrows_Avg,nchan_Avg,_,_=self.DicoData["flags_avg"].shape
        #     self.n4vis_Avg=self.n4vis_Avg_AllChan=nrows_Avg*nchan_Avg
        
        
        
        
        self.KernelMat_AllChan=None
        self.KernelMat_AllChan_Avg=None
        
        if shared_dict.exists(self.KernelSharedName):
            self.DicoKernelMats=shared_dict.attach(self.KernelSharedName)
            self.KernelMat_AllChan=self.DicoKernelMats.get("KernelMat_AllChan",None)
            if self.DoCompress:
                self.KernelMat_AllChan_Avg=self.DicoKernelMats.get("KernelMat_AllChan_Avg",None)
                Mask_Predict_Avg=self.DicoKernelMats.get("Mask_Predict_Avg",None)
                if self.KernelMat_AllChan_Avg is None:
                    stop
        else:
            self.DicoKernelMats=shared_dict.create(self.KernelSharedName)
            
        # self.DicoKernelMats=shared_dict.attach(self.KernelSharedName)
        # self.KernelMat_AllChan=self.DicoKernelMats.get("KernelMat_AllChan",None)
        # if self.DoCompress:
        #     self.KernelMat_AllChan_Avg=self.DicoKernelMats.get("KernelSharedNameAvg",None)

        if self.KernelMat_AllChan is not None or self.KernelMat_AllChan_Avg is not None:
            self.HasKernelMatrix=True
            if self.PolMode=="IFull":
                self.K_XX_AllChan=self.KernelMat_AllChan[0]
                self.K_YY_AllChan=self.KernelMat_AllChan[1]
                self.NJacobBlocks_X=2
                self.NJacobBlocks_Y=2
            elif self.PolMode=="Scalar":
                # n4vis=self.DicoData["data_flat"].size
                if self.KernelMat_AllChan is not None:
                    self.K_XX_AllChan=self.KernelMat_AllChan[0]
                    self.K_YY_AllChan=self.K_XX_AllChan
                    
                if self.DoCompress:
                    self.K_XX_AllChan_Avg=self.KernelMat_AllChan_Avg[0]
                    self.K_YY_AllChan_Avg=self.K_XX_AllChan_Avg
                    self.n4vis_Avg=self.n4vis_Avg_AllChan=self.NDirAvg * self.DicoData["NpBlBlocks"][0]

                # self.n4vis=n4vis
                self.NJacobBlocks_X=1
                self.NJacobBlocks_Y=1
            elif self.PolMode=="IDiag":
                # n4vis=self.DicoData["data_flat"].size
                self.K_XX_AllChan=self.KernelMat_AllChan[0]
                self.K_YY_AllChan=self.KernelMat_AllChan[1]
                # self.n4vis=n4vis
                self.NJacobBlocks_X=2
                self.NJacobBlocks_Y=1
            # self.Data=self.Data.reshape((nrows,nchan,self.NJacobBlocks,self.NJacobBlocks))
            # print("Kernel From shared")
            return
        else:
            # print("    COMPUTE KERNEL")
            pass

        T.timeit("stuff 2")
        # GiveArray(Name)
        nchan_AllChan=self.DicoData["freqs_full"].size
        n4vis_AllChan=nrows*nchan_AllChan
        self.n4vis_AllChan=n4vis_AllChan
            
        if self.PolMode=="IFull":
            #self.K_XX=np.zeros((NDir,n4vis/nchan,nchan),np.complex64)
            #self.K_YY=np.zeros((NDir,n4vis/nchan,nchan),np.complex64)
            self.DicoKernelMats["KernelMat_AllChan"]=np.zeros((2,NDir,n4vis_AllChan//nchan_AllChan,nchan_AllChan),dtype=self.CType)
            self.KernelMat_AllChan=self.DicoKernelMats["KernelMat_AllChan"]
            self.K_XX_AllChan=self.KernelMat_AllChan[0]
            self.K_YY_AllChan=self.KernelMat_AllChan[1]
            # KernelMatrix=NpShared.zeros(KernelSharedName,(n4vis,NDir,2),dtype=np.complex64)
            self.NJacobBlocks_X=2
            self.NJacobBlocks_Y=2
        elif self.PolMode=="Scalar":
            #n4vis=self.Data.size
            # KernelMatrix_XX=np.zeros((NDir,n4vis,nchan),np.complex64)
            # KernelMatrix=NpShared.zeros(KernelSharedName,(n4vis,NDir,1),dtype=np.complex64)
            self.DicoKernelMats["KernelMat_AllChan"]=np.zeros((1,NDir,n4vis_AllChan//nchan_AllChan,nchan_AllChan),dtype=self.CType)
            self.KernelMat_AllChan=self.DicoKernelMats["KernelMat_AllChan"]
            self.K_XX_AllChan=self.KernelMat_AllChan[0]
            self.K_YY_AllChan=self.K_XX_AllChan
            self.NJacobBlocks_X=1
            self.NJacobBlocks_Y=1


            if self.DoCompress:
                
                self.DicoKernelMats["KernelMat_AllChan_Avg"]=np.zeros((1,NDir,self.NDirAvg * self.DicoData["NpBlBlocks"][0],self.NChanSol),dtype=self.CType)
                self.KernelMat_AllChan_Avg=self.DicoKernelMats["KernelMat_AllChan_Avg"]
                self.K_XX_AllChan_Avg=self.KernelMat_AllChan_Avg[0]
                self.K_YY_AllChan_Avg=self.K_XX_AllChan_Avg
                self.n4vis_Avg=self.n4vis_Avg_AllChan=self.NDirAvg * self.DicoData["NpBlBlocks"][0]


            
            ##self.KernelMat_AllChan=NpShared.zeros(KernelSharedName,(1,NDir,n4vis_AllChan/nchan_AllChan,nchan_AllChan),dtype=self.CType)

            
        elif self.PolMode=="IDiag":
            self.DicoKernelMats["KernelMat_AllChan"]=np.zeros((2,NDir,n4vis_AllChan//nchan_AllChan,nchan_AllChan),dtype=self.CType)
            self.KernelMat_AllChan=self.DicoKernelMats["KernelMat_AllChan"]
            self.K_XX_AllChan=self.KernelMat_AllChan[0]
            self.K_YY_AllChan=self.KernelMat_AllChan[1]
            self.NJacobBlocks_X=2
            self.NJacobBlocks_Y=1
        T.timeit("stuff 3")
            
        #self.Data=self.Data.reshape((nrows,nchan,self.NJacobBlocks,self.NJacobBlocks))

        #self.K_XX=[]
        #self.K_YY=[]

        ApplyTimeJones=None
        #print(self.DicoData.keys())

        if shared_dict.exists("DicoPreApplyJones_Chunk%i"%self.iTimeChunk):
            DicoPreApplyJones=shared_dict.attach("DicoPreApplyJones_Chunk%i"%self.iTimeChunk)
            if len(DicoPreApplyJones)>0:
                ApplyTimeJones=DicoPreApplyJones

        #import gc
        #gc.enable()
        # gc.set_debug(gc.DEBUG_LEAK)


        # ##############################################
        # from SkyModel.Sky import ClassSM
        # SM=ClassSM.ClassSM("ModelRandom00.txt.npy")
        # SM.Type="Catalog"
        # rac,decc=self.DATA["rac_decc"]
        # SM.Calc_LM(rac,decc)
        # # SM.SourceCat=SM.SourceCat[self.SM.MapClusterCatOrigToCut]
        # # SM.ClusterCat=SM.ClusterCat[self.SM.MapClusterCatOrigToCut]
        # SM.NDir=SM.SourceCat.shape[0]
        # SM.Dirs=SM.Dirs
        
        # # Namemm="SM_Compress_%i.npy"%self.iAnt
        # # np.save(Namemm,self.SM_Compress.SourceCat)
        # # SM=ClassSM.ClassSM(Namemm)
        # # #SM.SourceCat.I*=10
        # # #SM.SourceCat.Sref*=10
        # # SM.Calc_LM(rac,decc)
        
        # self.KernelMat1=np.zeros((1,NDir,n4vis//nchan,nchan),dtype=self.CType)
        # self.K1_XX=self.KernelMat1[0]
        # self.K1_YY=self.K1_XX

        iAntSel=23
        
        # if self.iAnt==iAntSel:
        #     import pylab
        #     pylab.figure(0)
        #     pylab.clf()
            
        # pylab.figure(1)
        # pylab.clf()
        # pylab.figure(2)
        # pylab.clf()

        


        for iDir in range(NDir):
            if self.SM.Type=="Column":
                K=self.DicoData["data_predict"]
            elif False:#"PredictArray" in self.DATA.keys():
                Kall=self.DATA["PredictArray"][iDir]
                ind0=self.DicoData["indOrig"]
                ind1=self.DicoData["indOrig1"]
                D0=Kall[ind0]
                D1=Kall[ind1].conj()
                c1=D1[:,:,1].copy()
                c2=D1[:,:,2].copy()
                D1[:,:,1]=c2
                D1[:,:,2]=c1
                K=np.concatenate([D0,D1],axis=0)
                # uvw1=self.DATA["uvw"]
                # uvw1=np.concatenate([uvw1[ind0],-uvw1[ind1]],axis=0)
                # uvw0=self.DicoData["uvw"]
                # K0=self.PM.predictKernelPolCluster(self.DicoData,self.SM,iDirection=iDir,ApplyTimeJones=ApplyTimeJones)
            elif self.SM.Type=="Hybrid":
                ThisSM,ThisIDir=self.SM.give_SM_iDir(iDir)
                #print("K: ",iDir,ThisSM.Type,ThisIDir)
                K=self.PM.predictKernelPolCluster(self.DicoData,ThisSM,iDirection=ThisIDir,ApplyTimeJones=ApplyTimeJones)#,iDirJones=ThisIDir)
            else:
                K=self.PM.predictKernelPolCluster(self.DicoData,self.SM,iDirection=iDir,ApplyTimeJones=ApplyTimeJones)
                
            if np.count_nonzero(np.isnan(K))>1:
                print(self.iAnt,iDir,K)
                stop

            #K=self.PM.predictKernelPolCluster(self.DicoData,self.SM,iDirection=iDir)#,ApplyTimeJones=ApplyTimeJones)
            #K*=-1
            T.timeit("Calc K0")

            
                #gc.collect()
                #print(gc.garbage)


            # if (iDir==31)&(self.iAnt==51):
            #     ifile=0
            #     while True:
            #         fname="png/Kernel.%5.5i.npy"%ifile
            #         if not(os.path.isfile(fname)) :
            #             np.save(fname,K)
            #             break
            #         ifile+=1

            K_XX=K[:,:,0]
            K_YY=K[:,:,3]
            if self.PolMode=="Scalar":
                K_XX=(K_XX+K_YY)/2.
                K_YY=K_XX

            self.K_XX_AllChan[iDir,:,:]=K_XX
            self.K_YY_AllChan[iDir,:,:]=K_YY
            #self.K_XX.append(K_XX)
            #self.K_YY.append(K_YY)


        #     ######################
        #     if self.iAnt==iAntSel:
        #         # K1=self.PM.predictKernelPolCluster(self.DicoData,SM,iDirection=iDir)#,ApplyTimeJones=ApplyTimeJones)
        #         # # K1=self.PM.predictKernelPolCluster(self.DicoData,self.SM,iDirection=iDir,ApplyTimeJones=ApplyTimeJones,ForceNoDecorr=True)
        #         # K1[K==0]=0
        #         # K1_XX=K1[:,:,0]
        #         # K1_YY=K1[:,:,3]
        #         # if self.PolMode=="Scalar":
        #         #     K1_XX=(K1_XX+K1_YY)/2.
        #         #     K1_YY=K1_XX
        #         # self.K1_XX[iDir,:,:]=K1_XX
        #         # self.K1_YY[iDir,:,:]=K1_YY
                
        #         A0=self.DicoData["A0"]
        #         A1=self.DicoData["A1"]
        #         ind=np.arange(K.shape[0])#np.where((A0==0)&(A1==26))[0]
        #         d1=K[ind,0,0] 
        #         # ind=np.arange(K1.shape[0])#np.where((A0==0)&(A1==26))[0]
        #         # d0=K1[ind,0,0]
        #         op0=np.abs
        #         #op0=np.real
                
        #         op1=np.imag
        #         pylab.figure(0)
        #         #print(iDir)
        #         #pylab.subplot(1,NDir,iDir+1)
        #         #pylab.plot(op0(d0))
        #         color="gray"
        #         if iDir==26:
        #             color="red"
        #             GGG=op0(d1)
        #         pylab.plot(op0(d1),
        #                    alpha=0.1,
        #                    color=color)
        #         # # pylab.plot(op0(d0)-op0(d1))
        #         # # pylab.plot(op0(d1)/op0(d0))
        #         # #pylab.ylim(-100,100)
        #         # pylab.draw()
        #         # pylab.show(block=False)
                
        #         # # op1=np.angle
        #         # # pylab.figure(1)
        #         # # pylab.subplot(1,NDir,iDir+1)
        #         # # pylab.plot(op1(d0))
        #         # # pylab.plot(op1(d1))
        #         # # pylab.plot(op1(d0)-op1(d1))
        #         # # #pylab.plot(op1(d1)-op1(d0))
        #         # # pylab.ylim(-np.pi,np.pi)
                
        #         # # pylab.subplot(2,1,2)
        #         # # #pylab.plot(op1(d0))
        #         # # pylab.plot(op1(d1*d0.conj()))#,ls="--")
        #         # # #pylab.plot(op1(d0*d1.conj()),ls="--")
        #         # # #pylab.ylim(-1,1)
        #         # pylab.draw()
        #         # pylab.show(block=False)
        #         # pylab.pause(0.1)

        if self.DoCompress:
            T.reinit()
            k=self.AverageMachine.AverageKernelMatrix(self.DicoData,self.K_XX_AllChan)
            T.timeit("AverageKernelMatrix")
            
            self.K_XX_AllChan_Avg[:,:,:]=k[:,:,:]
            self.K_YY_AllChan_Avg[:]=self.K_XX_AllChan_Avg[:]
            #print("iAnt",self.iAnt)
            Mask_Predict_Avg=(self.K_XX_AllChan==0)

            LMask_Predict_Avg=[]
            for iChanSol in range(self.NChanSol):
                LMask_Predict_Avg.append(np.any(self.K_XX_AllChan[:,:,(self.VisToSolsChanMapping==iChanSol)]==0,axis=0))
            Mask_Predict_Avg=np.array(LMask_Predict_Avg)
            self.DicoKernelMats["Mask_Predict_Avg"]=Mask_Predict_Avg
            self.DicoKernelMats.delete_item("KernelMat_AllChan")
            
        # if self.iAnt==iAntSel:
        #     pylab.plot(GGG,
        #                color="red")
        #     pylab.draw()
        #     pylab.show(block=False)
        #     pylab.pause(0.1)
            
                

                
        # #     del(K1,K1_XX,K1_YY)
        # #     del(K,K_XX,K_YY)
        # stop


        # 
        #stop
        #gc.collect()

        #print("DicoToShared %s"%self.SharedDataDicoName)
        
        self.HasKernelMatrix=True
        T.timeit("stuff 4")


    def AntennaDicoData_to_shm(self):
        if shared_dict.exists(self.SharedDataDicoName):
            return
        
        if self.DoCompress:
            T=ClassTimeIt.ClassTimeIt("  AntennaDicoData_to_shm")
            T.disable()
            T.reinit()
            Mask_Predict_Avg=self.DicoKernelMats["Mask_Predict_Avg"]

            ch0,ch1=self.DicoData["ch0ch1"]
            iChanSol0=self.VisToSolsChanMapping[ch0]
            iChanSol1=self.VisToSolsChanMapping[ch1-1]
            if iChanSol1!=iChanSol0: stop
            iChanSol=iChanSol1
            self.AverageMachine.AverageDataVector(self.DicoData,
                                                  Mask=Mask_Predict_Avg[iChanSol],
                                                  Stop=(self.iAnt==13))
            T.timeit("AverageDataVector")
            del(self.DicoData['data'],self.DicoData['data_flat'])
            
            # if True:#self.iAnt==13:
                
            #     op0=np.abs
                
            #     #k1=self.AverageMachine.AverageKernelMatrix(self.DicoData,self.K1_XX)
                
            #     # pylab.figure(2)
            #     # for iDir in range(NDir):
            #     #     d0=k[iDir,:,0]
            #     #     d1=k1[iDir,:,0]
            #     #     pylab.subplot(1,NDir,iDir+1)
            #     #     pylab.plot(op0(d0))
            #     #     pylab.plot(op0(d1))
            #     #     pylab.plot(op0(d0)-op0(d1))
            #     #     pylab.legend(("k0","k1","rk"))
            #     #     #pylab.plot(op0(d0)-op0(d1))
            #     #     #pylab.plot(op0(d1)/op0(d0))
            #     # #pylab.ylim(-100,100)
            #     # pylab.draw()
            #     # pylab.show(block=False)
                
                
            #     f=np.any(k[:,:,0]==0,axis=0)
            #     import pylab
            #     pylab.figure(3)
            #     d0=np.sum(k[:,:,0],axis=0)
            #     #d1=np.sum(k1[:,:,0],axis=0)
            #     dd=self.DicoData["data_flat_avg"][0]
            #     # d0.flat[f]=0
            #     # d1.flat[f]=0
            #     pylab.clf()
            #     pylab.plot(op0(d0))
            #     #pylab.plot(op0(d1))
            #     pylab.plot(op0(dd))
            #     pylab.plot(op0(dd)-op0(d0))
            #     #pylab.legend(("k0","k1","d","r"))
            #     pylab.legend(("k0","d","r"))
            #     pylab.title(self.iAnt)
            #     pylab.draw()
            #     pylab.show(block=False)
            #     pylab.pause(1)
        self.DicoData=shared_dict.dict_to_shm(self.SharedDataDicoName,self.DicoData)


        
    def SelectChannelKernelMat(self):

        if self.DoCompress:
            flags_key="flags_avg"
            flags_flat_key="flags_flat_avg"
            ch0,ch1=self.DicoData["ch0ch1"]
            iChanSol0=self.VisToSolsChanMapping[ch0]
            iChanSol1=self.VisToSolsChanMapping[ch1-1]
            if iChanSol1!=iChanSol0: stop
            self.K_XX=self.K_XX_AllChan_Avg[:,:,iChanSol0:iChanSol1+1]
            self.K_YY=self.K_YY_AllChan_Avg[:,:,iChanSol0:iChanSol1+1]
            if self.DoMergeStations: return
        else:
            self.K_XX=self.K_XX_AllChan[:,:,self.ch0:self.ch1]
            self.K_YY=self.K_YY_AllChan[:,:,self.ch0:self.ch1]
            flags_key="flags"
            flags_flat_key="flags_flat"
            
        flags_orig_key="flags"
        flags_orig_flat_key="flags_flat"
        
        DicoData=self.DicoData
        # print(DicoData["flags"].shape,DicoData["data"].shape,DicoData["flags_flat"].shape,DicoData["data_flat"].shape,self.K_XX.shape)
        # if self.DoCompress:
        #     print(DicoData["flags_avg"].shape,DicoData["data_avg"].shape,DicoData["flags_flat_avg"].shape,DicoData["data_flat_avg"].shape)

        # No-Comp:
        # anal
        # (840, 4, 1, 1) (840, 4, 1, 1) (1, 3360) (1, 3360) (9, 840, 4)
        # fft:
        # (840, 4, 1, 1) (840, 4, 1, 1) (1, 3360) (1, 3360) (9, 840, 4)
        # Compress
        # anal    :
        #(840, 4, 1, 1) (840, 4, 1, 1) (1, 3360) (1, 3360) (9, 189, 1)
        #(189, 1, 1, 1) (189, 1, 1, 1) (1, 189) (1, 189)
        # fft:
        #(840, 4, 1, 1) (840, 4, 1, 1) (1, 3360) (1, 3360) (9, 189, 1)
        #(189, 1, 1, 1) (189, 1, 1, 1) (1, 189) (1, 189)

        NDir=self.SM.NDir
        _,nr,nch=self.K_XX.shape
        for iDir in range(NDir):
            if self.SM.Type=="Hybrid":
                SM,ThisIDir=self.SM.give_SM_iDir(iDir)
            else:
                SM=self.SM
                ThisIDir=iDir
                
            indSources=np.where(SM.SourceCat.Cluster==ThisIDir)[0]
            Type=np.unique(SM.SourceCat.Type[indSources])
            if not 0 in Type:
                #print("!!!!!!!!::::::::::::::")
                continue
            
            K=self.K_XX[iDir,:,:]

            indRow,indChan=np.where(K==0)

            # print("")
            # print("non-zeroflag:",indRow.size)
            # print("")
            # f=self.DicoData[flags_key]
            # f0=np.count_nonzero(f==1)/f.size
            self.DicoData[flags_key][indRow,indChan,:]=1
            # f1=np.count_nonzero(f==1)/f.size
            # print("[%i] %f -> %f"%(self.iAnt,f0,f1))
            #self.DicoData[flags_key].reshape((NDir*nr,nch,1))[indRow,indChan,:]=1
            #stop
            #self.DicoData[flags_orig_key][indRow,indChan,:]=1
            

            
        DicoData=self.DicoData
        nr,nch = K.shape
        # print(self.iAnt,DicoData["data"].shape,DicoData[flags_key].shape,self.NJacobBlocks_X,nr,nch,self.NJacobBlocks_Y)
        flags_flat = np.rollaxis(DicoData[flags_key],2).reshape(self.NJacobBlocks_X,nr*nch*self.NJacobBlocks_Y)
        DicoData[flags_flat_key][flags_flat]=1
        #nr,nch = self.K_YY_AllChan[0,:,self.ch0:self.ch1]#K.shape
        #flags_flat=np.rollaxis(DicoData[flags_orig_key],2).reshape(self.NJacobBlocks_X,nr*nch*self.NJacobBlocks_Y)
        #DicoData[flags_orig_flat_key][flags_flat]=1


        self.DataAllFlagged=False
        NP,_=DicoData[flags_flat_key].shape
        for ipol in range(NP):
            f=(DicoData[flags_flat_key][ipol]==0)
            ind=np.where(f)[0]
            if ind.size==0: 
                self.DataAllFlagged=True
                continue
            fracFlagged=ind.size/float(f.size)
            if fracFlagged<self.maxFlagFrac:#ind.size==0:
                # print(self.iAnt,"DataAllFlagged",fracFlagged)
                self.DataAllFlagged=True

        # if self.DataAllFlagged:
        #     stop
        # # print("ok!!!")
        # # stop

        # print("SelectChannelKernelMat",np.count_nonzero(DicoData["flags_flat"]),np.count_nonzero(DicoData["flags"]))



    def GiveData(self,DATA,iAnt,rms=0.):
        
        #DicoData=NpShared.SharedToDico(self.SharedDataDicoName)

        DicoData={}
        if shared_dict.exists(self.SharedDataDicoName):
            DicoData=shared_dict.attach(self.SharedDataDicoName)
            #print("%s exists"%self.SharedDataDicoName)
        else:
            #print("%s not exists"%self.SharedDataDicoName)
            pass
            
        if len(DicoData)==0:
            #print("     COMPUTE DATA")
            DicoData={}
            ind0=np.where(DATA['A0']==iAnt)[0]
            ind1=np.where(DATA['A1']==iAnt)[0]
            self.ZeroSizedData=False
            if ind0.size==0 and ind1.size==0:
                self.ZeroSizedData=True
            DicoData["A0"] = np.concatenate([DATA['A0'][ind0], DATA['A1'][ind1]])
            DicoData["A1"] = np.concatenate([DATA['A1'][ind0], DATA['A0'][ind1]])
            NRow,Nch_all,NPol=DATA['data'].shape
            #ChanIndex=np.ones((NRow,1,1))*np.int8(np.arange(Nch_all).reshape((1,-1,1)))
            DicoData["ch0ch1"]=(self.ch0,self.ch1)
            D0=DATA['data'][ind0,self.ch0:self.ch1]
            D1=DATA['data'][ind1,self.ch0:self.ch1].conj()

            #ChanIndex0=ChanIndex[ind0,:]
            #ChanIndex1=ChanIndex[ind1,:]
            
            A0,A1=DicoData["A0"],DicoData["A1"]
            A0A1=sorted(list(set([(A0[i],A1[i]) for i in range(DicoData["A0"].size)])))
            NpBlBlocks=len(A0A1)
            DicoData["NpBlBlocks"]=np.array([NpBlBlocks])


            
            c1=D1[:,:,1].copy()
            c2=D1[:,:,2].copy()
            D1[:,:,1]=c2
            D1[:,:,2]=c1

            
            
            DicoData["data"] = np.concatenate([D0, D1])
            #DicoData["ChanIndex"] = np.concatenate([ChanIndex0, ChanIndex1])
            if self.SM.Type=="Column":

                # print("=================")
                # print(self.DicoMainDataName)
                # print(self.DicoColPredictDataName)
                # os.system("ls /dev/shm/ddf.*/%s*"%self.DicoMainDataName)
                # os.system("ls /dev/shm/ddf.*/%s*"%self.DicoColPredictDataName)
                # print(list(self.DATA_ColumnPredict.keys()))
                # print(self.DATA_ColumnPredict.path)
                # print("=================")
                
                data_predict=self.DATA_ColumnPredict['data'].copy()
                #print("OKKK")
                D0=data_predict[ind0,:]
                D1=data_predict[ind1,:].conj()
                c1=D1[:,:,1].copy()
                c2=D1[:,:,2].copy()
                D1[:,:,1]=c2
                D1[:,:,2]=c1
                DicoData["data_predict"] = np.concatenate([D0, D1])


            DicoData["indOrig"] = ind0
            DicoData["indThis"]=np.arange(DicoData["indOrig"].size)
            DicoData["indOrig1"] = ind1
            DicoData["uvw"]  = np.concatenate([DATA['uvw'][ind0], -DATA['uvw'][ind1]])
            if "UVW_dt" in DATA.keys():
                DicoData["UVW_dt"]  = np.concatenate([DATA["UVW_dt"][ind0], -DATA["UVW_dt"][ind1]])

            if "W" in DATA.keys():
                DicoData["W"] = np.concatenate([DATA['W'][ind0,self.ch0:self.ch1], DATA['W'][ind1,self.ch0:self.ch1]])

            # DicoData["IndexTimesThisChunk"]=np.concatenate([DATA["IndexTimesThisChunk"][ind0], DATA["IndexTimesThisChunk"][ind1]]) 
            # DicoData["UVW_RefAnt"]=DATA["UVW_RefAnt"][it0:it1]

            if "Kp" in DATA.keys():
                 DicoData["Kp"]=DATA["Kp"]

            D0=DATA['flags'][ind0,self.ch0:self.ch1]
            D1=DATA['flags'][ind1,self.ch0:self.ch1].conj()
            c1=D1[:,:,1].copy()
            c2=D1[:,:,2].copy()
            D1[:,:,1]=c2
            D1[:,:,2]=c1
            DicoData["flags"] = np.concatenate([D0, D1])

            D0=DATA['flags'][ind0,:]
            D1=DATA['flags'][ind1,:].conj()
            c1=D1[:,:,1].copy()
            c2=D1[:,:,2].copy()
            D1[:,:,1]=c2
            D1[:,:,2]=c1
            DicoData["flags_allFreqs"] = np.concatenate([D0, D1])
            
            if self.SM.Type=="Image" or self.SM.Type=="Hybrid":
                #DicoData["flags_image"]=DicoData["flags"].copy()
                nr,_,_=DicoData["data"].shape
                _,nch,_=DATA['data'].shape
                DicoData["flags_image"]=np.zeros((nr,nch,4),np.bool8)
                #DicoData["flags_image"].fill(0)

            nr,nch,_=DicoData["data"].shape

            if self.GD["VisData"]["FreePredictGainColName"]!=None or self.GD["VisData"]["FreePredictColName"]!=None:
                DicoData["data_allpol"]=DicoData["data"].reshape((nr,nch,2,2))
                
            if self.PolMode=="Scalar":
                d=(DicoData["data"][:,:,0]+DicoData["data"][:,:,-1])/2
                DicoData["data"] = d.reshape((nr,nch,1))

                # if self.SM.Type=="Column":
                #     d=(DicoData["data_predict"][:,:,0]+DicoData["data_predict"][:,:,-1])/2
                #     DicoData["data_predict"] = d.reshape((nr,nch,1))


                f=(DicoData["flags"][:,:,0]|DicoData["flags"][:,:,-1])
                DicoData["flags"] = f.reshape((nr,nch,1))

                f=(DicoData["flags_allFreqs"][:,:,0]|DicoData["flags_allFreqs"][:,:,-1])
                DicoData["flags_allFreqs"] = f.reshape((nr,Nch_all,1))

                
                if self.GD["Compression"]["CompressionMode"]:
                    DicoData["A0_freq"]=(DicoData["A0"].reshape((nr,1,1)) * np.ones((1,nch,1)))
                    DicoData["A1_freq"]=(DicoData["A1"].reshape((nr,1,1)) * np.ones((1,nch,1)))

            elif self.PolMode=="IDiag":
                d=DicoData["data"][:,:,0::3]
                DicoData["data"] = d.copy().reshape((nr,nch,2))

                # if self.SM.Type=="Column":
                #     d=DicoData["data_predict"][:,:,0::3]
                #     DicoData["data_predict"] = d.copy().reshape((nr,nch,2))

                f=DicoData["flags"][:,:,0::3]
                DicoData["flags"] = f.copy().reshape((nr,nch,2))

                f=DicoData["flags_allFreqs"][:,:,0::3]
                DicoData["flags_allFreqs"] = f.reshape((nr,Nch_all,2))


            DicoData["freqs"]   = DATA['freqs'][self.ch0:self.ch1]
            DicoData["dfreqs"]   = DATA['dfreqs'][self.ch0:self.ch1]
            DicoData["times"] = np.concatenate([DATA['times'][ind0], DATA['times'][ind1]])
            DicoData["infos"] = DATA['infos']

            # nr,nch,_=DicoData["data"].shape

            FlagsShape=DicoData["flags"].shape
            FlagsSize=DicoData["flags"].size
            
            DicoData["flags"]=DicoData["flags"].reshape(nr,nch,self.NJacobBlocks_X,self.NJacobBlocks_Y)
            DicoData["flags_allFreqs"]=DicoData["flags_allFreqs"].reshape(nr,Nch_all,self.NJacobBlocks_X,self.NJacobBlocks_Y)
            DicoData["data"]=DicoData["data"].reshape(nr,nch,self.NJacobBlocks_X,self.NJacobBlocks_Y)
            # if self.SM.Type=="Column":
            #     DicoData["data_predict"]=DicoData["data_predict"].reshape(nr,nch,self.NJacobBlocks_X,self.NJacobBlocks_Y)

            DicoData["flags_flat"]=np.rollaxis(DicoData["flags"],2).reshape(self.NJacobBlocks_X,nr*nch*self.NJacobBlocks_Y)
            DicoData["flags_allFreqs_flat"]=np.rollaxis(DicoData["flags_allFreqs"],2).reshape(self.NJacobBlocks_X,nr*Nch_all*self.NJacobBlocks_Y)
            DicoData["data_flat"]=np.rollaxis(DicoData["data"],2).reshape(self.NJacobBlocks_X,nr*nch*self.NJacobBlocks_Y)
            #DicoData["ChanIndex_flat"]=np.rollaxis(DicoData["ChanIndex"],2).reshape(self.NJacobBlocks_X,nr*Nch_all*self.NJacobBlocks_Y)


            if self.GD["Compression"]["CompressionMode"] is not None:
                DicoData["A0_freq_flat"]=np.rollaxis(DicoData["A0_freq"],2).reshape(self.NJacobBlocks_X,nr*nch*self.NJacobBlocks_Y)
                DicoData["A1_freq_flat"]=np.rollaxis(DicoData["A1_freq"],2).reshape(self.NJacobBlocks_X,nr*nch*self.NJacobBlocks_Y)


            # DicoData["data_predict_flat"]=np.rollaxis(DicoData["data_predict"],2).reshape(self.NJacobBlocks_X,nr*nch*self.NJacobBlocks_Y)
            # ###################
            # NJacobBlocks_X=2
            # NJacobBlocks_Y=2
            # F0=np.zeros((nr,nch,NJacobBlocks_X,NJacobBlocks_Y))
            # FlagsShape=F0.shape
            # FlagsSize=F0.size
            # F0=np.arange(FlagsSize).reshape(FlagsShape)
            # F0Flat=np.rollaxis(F0,2).reshape(NJacobBlocks_X,nr*nch*NJacobBlocks_Y)
            # F1=np.rollaxis(F0Flat.reshape(NJacobBlocks_X,nr,nch,NJacobBlocks_Y),0,3).reshape(FlagsShape)
            # print(np.count_nonzero((F0-F1).ravel()))
            # stop
            # ###################


            # del(DicoData["data"])

            C0=("W" in DicoData.keys())
            C_KAFCA=(self.GD["Solvers"]["SolverType"]=="KAFCA")
            C_CohJones=(self.GD["Solvers"]["SolverType"]=="CohJones")
            C_CohJones_with_weights = C_CohJones and (self.GD["Weighting"]["WeightInCol"] is not None) and (self.GD["Weighting"]["WeightInCol"]!="")
            if C_CohJones_with_weights:
                rms=1.
            C_HasData = (not self.ZeroSizedData)
            self.WillUseWeights=C_HasData and (C_KAFCA or C_CohJones_with_weights)
            

            # if not self.WillUseWeights:
            #     print("Not using weihts",rms)
                
            if self.WillUseWeights:
                #print("Using weihts",rms)
                DicoData["rms"]=np.array([rms],np.float32)
                u,v,w=DicoData["uvw"].T
                if self.ResolutionRad!=None:
                    freqs=DicoData["freqs"]
                    wave=np.mean(299792456./freqs)
                    d=np.sqrt((u/wave)**2+(v/wave)**2)
                    FWHMFact=2.*np.sqrt(2.*np.log(2.))
                    sig=self.ResolutionRad/FWHMFact
                    V=(1./np.exp(-d**2*np.pi*sig**2))**2
                    V=V.reshape((V.size,1,1))*np.ones((1,freqs.size,self.npolData))
                else:
                    V=np.ones((u.size,freqs.size,self.npolData),np.float32)
                    
                if "W" in DicoData.keys():
                    W=DicoData["W"]**2
                    W_nrows,W_nch=W.shape
                    W[W==0]=1.e-6
                    V=V/W.reshape((W_nrows,W_nch,1))

                R=rms**2*V
                #print(R.min(),R.max())
                
                Rinv=1./R
                #print(Rinv.min(),Rinv.max())
                
                Weights=W.reshape((W_nrows,W_nch,1))*np.ones((1,1,self.npolData))

                # DicoData["Rinv"]=Rinv
                # DicoData["R"]=R
                # DicoData["Weights"]=Weights
                
                self.R_flat=np.rollaxis(R,2).reshape(self.NJacobBlocks_X,nr*nch*self.NJacobBlocks_Y)
                self.Rinv_flat=np.rollaxis(Rinv,2).reshape(self.NJacobBlocks_X,nr*nch*self.NJacobBlocks_Y)
                self.Weights_flat=np.rollaxis(Weights,2).reshape(self.NJacobBlocks_X,nr*nch*self.NJacobBlocks_Y)

                self.R_flat=np.require(self.R_flat,dtype=self.CType)
                self.Rinv_flat=np.require(self.Rinv_flat,dtype=self.CType)
                
                Rmin=np.min(R)
                #Rmax=np.max(R)
                Flag=(self.R_flat>1e3*Rmin)
                DicoData["flags_flat"][Flag]=1
                DicoData["Rinv_flat"]=self.Rinv_flat
                DicoData["R_flat"]=self.R_flat
                DicoData["Weights_flat"]=self.Weights_flat
                

            self.DataAllFlagged=False
            NP,_=DicoData["flags_flat"].shape
            for ipol in range(NP):
                f=(DicoData["flags_flat"][ipol]==0)
                ind=np.where(f)[0]
                
                if ind.size==0: 
                    self.DataAllFlagged=True
                    continue

                fracFlagged=ind.size/float(f.size)
                if fracFlagged<self.maxFlagFrac:#ind.size==0:
                    self.DataAllFlagged=True

            # if self.DataAllFlagged:
            #     stop
            # # DicoData=NpShared.DicoToShared(self.SharedDataDicoName,DicoData)
            # # self.SharedDicoDescriptors["SharedAntennaVis"]=NpShared.SharedDicoDescriptor(self.SharedDataDicoName,DicoData)


            if shared_dict.exists("DicoPreApplyJones_Chunk%i"%self.iTimeChunk):
                DicoPreApplyJones=shared_dict.attach("DicoPreApplyJones_Chunk%i"%self.iTimeChunk)
                DicoJonesMatrices={}
                ind0=DicoData["indOrig"]
                ind1=DicoData["indOrig1"]
                MapTimes=DATA["Map_VisToJones_Time"]
                MapTimesSel=np.concatenate([MapTimes[ind0], MapTimes[ind1]])
                DicoData["Map_VisToJones_Time"]=MapTimesSel
                DicoData["Map_VisToJones_Freq"]=DATA["Map_VisToJones_Freq"]
                
            DicoData["freqs_full"]   = self.DATA['freqs']
            DicoData["dfreqs_full"]   = self.DATA['dfreqs']

            

        # Antenna based Dico is already computed in a shared 
        else:
            
            self.ZeroSizedData=False
            if DicoData["A0"].size==0 and DicoData["A1"].size==0:
                self.ZeroSizedData=True

            C0=("W" in DicoData.keys())
            C_KAFCA=(self.GD["Solvers"]["SolverType"]=="KAFCA")
            C_CohJones=(self.GD["Solvers"]["SolverType"]=="CohJones")
            C_CohJones_with_weights = C_CohJones and (self.GD["Weighting"]["WeightInCol"] is not None) and (self.GD["Weighting"]["WeightInCol"]!="")
            if C_CohJones_with_weights:
                rms=1.
            C_HasData = (not self.ZeroSizedData)
            self.WillUseWeights=C_HasData and (C_KAFCA or C_CohJones_with_weights)

                
            if self.WillUseWeights:
                self.Rinv_flat=DicoData["Rinv_flat"]
                self.R_flat=DicoData["R_flat"]
                self.Weights_flat=DicoData["Weights_flat"]
            
            #print("DATA From shared")
            #print(np.max(DicoData["A0"]))
            #np.save("testA0",DicoData["A0"])
            #DicoData["A0"]=np.load("testA0.npy")
            #DicoData=NpShared.SharedToDico(self.SharedDataDicoName)
            #print(np.max(DicoData["A0"]))
            #print
            #stop

        


            
            #print(DATA["Map_VisToJones_Time"].max())
            #stop

        self.DoTikhonov=False
        #self.GD["CohJones"]["LambdaTk"]=0
        if (self.GD["CohJones"]["LambdaTk"]!=0)&(self.GD["Solvers"]["SolverType"]=="CohJones"):
            self.DoTikhonov=True
            self.LambdaTk=self.GD["CohJones"]["LambdaTk"]
            self.DicoTikhonov=shared_dict.attach("DicoTikhonov")
            self.Linv=self.DicoTikhonov["Linv"]
            self.X0=self.DicoTikhonov["X0"]
            

        # DicoData["A0"] = np.concatenate([DATA['A0'][ind0]])
        # DicoData["A1"] = np.concatenate([DATA['A1'][ind0]])
        # D0=DATA['data'][ind0]
        # DicoData["data"] = np.concatenate([D0])
        # DicoData["uvw"]  = np.concatenate([DATA['uvw'][ind0]])
        # DicoData["flags"] = np.concatenate([DATA['flags'][ind0]])
        # DicoData["freqs"]   = DATA['freqs']

        self.DataAllFlagged=False
        try:
            NP,_=DicoData["flags_flat"].shape
        except:
            print(DicoData.keys())
            stop
            
        for ipol in range(NP):
            f=(DicoData["flags_flat"][ipol]==0)
            ind=np.where(f)[0]
            if ind.size==0: 
                self.DataAllFlagged=True
                continue
            fracFlagged=ind.size/float(f.size)
            if fracFlagged<self.maxFlagFrac:#ind.size==0:
                self.DataAllFlagged=True


        # if self.DoCompress:
        #     self.AverageMachine.AverageDataVector(DicoData)

        return DicoData

###########################################
###########################################
###########################################


def testPredict():
    import pylab
    VS=ClassVisServer.ClassVisServer("../TEST/0000.MS/")

    MS=VS.MS
    SM=ClassSM.ClassSM("../TEST/ModelRandom00.txt.npy")
    SM.Calc_LM(MS.rac,MS.decc)




    nd=SM.NDir
    npol=4
    na=MS.na
    Gains=np.zeros((na,nd,npol),dtype=np.complex64)
    Gains[:,:,0]=1
    Gains[:,:,-1]=1
    #Gains+=np.random.randn(*Gains.shape)*0.5+1j*np.random.randn(*Gains.shape)
    Gains=np.random.randn(*Gains.shape)+1j*np.random.randn(*Gains.shape)
    #Gains[:,1,:]=0
    #Gains[:,2,:]=0
    #g=np.random.randn(*(Gains[:,:,0].shape))+1j*np.random.randn(*(Gains[:,:,0].shape))
    #g=g.reshape((na,nd,1))
    #Gains*=g

    DATA=VS.GiveNextVis(0,50)

    # Apply Jones
    PM=ClassPredict(Precision="S",BeamAtFacet=(self.GD["Beam"]["BeamAt"].lower() == "facet"))
    DATA["data"]=PM.predictKernelPolCluster(DATA,SM,ApplyJones=Gains)
    
    ############################
    PolMode="IFull"#"Scalar"
    iAnt=10
    JM=ClassJacobianAntenna(SM,iAnt,PolMode=PolMode)
    JM.setDATA(DATA)
    JM.CalcKernelMatrix()
    if PolMode=="Scalar":
        Gains=Gains[:,:,0].reshape((na,nd,1))

    Jacob= JM.CalcJacobianAntenna(Gains)

    y=JM.GiveDataVec()
    
#    Gain=JM.ThisGain[:,1,:]
    predict=JM.J_x(Gains[iAnt])

    pylab.figure(1)
    pylab.clf()
    pylab.subplot(2,1,1)
    pylab.plot(predict.real)
    pylab.plot(y.real)
    pylab.plot((predict-y).real)
    pylab.subplot(2,1,2)
    pylab.plot(predict.imag)
    pylab.plot(y.imag)
    pylab.plot((predict-y).imag)
    pylab.draw()
    pylab.show(False)

    pylab.figure(2)
    pylab.clf()
    pylab.subplot(1,2,1)
    pylab.imshow(np.abs(JM.JHJ),interpolation="nearest")
    pylab.subplot(1,2,2)
    pylab.imshow(np.abs(JM.JHJinv),interpolation="nearest")
    pylab.draw()
    pylab.show(False)

    pylab.figure(3)
    pylab.clf()
    pylab.imshow(np.abs(JM.Jacob)[0:20],interpolation="nearest")
    pylab.draw()
    pylab.show(False)

    stop    
    
