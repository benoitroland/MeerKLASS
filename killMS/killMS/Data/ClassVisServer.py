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
from pyrap.tables import table
from DDFacet.Other import logger
global log
log=logger.getLogger("ClassVisServer")
# import MyPickle
#from killMS.Array import NpShared
from killMS.Other import ClassTimeIt
from killMS.Other import ModColor
from killMS.Array import ModLinAlg
#slogger.setSilent(["NpShared"])
#from Sky.PredictGaussPoints_NumExpr3 import ClassPredictParallel as ClassPredict 
#from Sky.PredictGaussPoints_NumExpr3 import ClassPredict as ClassPredict 
from . import ClassWeighting
from killMS.Other import reformat
import os
from killMS.Other.ModChanEquidistant import IsChanEquidistant
from killMS.Data import ClassReCluster
#import MergeJones
from killMS.Data import ClassJonesDomains
#from DDFacet.Imager import ClassWeighting as ClassWeightingDDF
#from DDFacet.Other.PrintList import ListToStr
import DDFacet.Other.PrintList
from killMS.Other import ClassGiveSolsFile
from DDFacet.Other import AsyncProcessPool as APP

from DDFacet.Array import shared_dict
#from DDFacet.Other.AsyncProcessPool import APP
from collections import OrderedDict
import time

def AngDist(ra0,ra1,dec0,dec1):
    AC=np.arccos
    C=np.cos
    S=np.sin
    D=S(dec0)*S(dec1)+C(dec0)*C(dec1)*C(ra0-ra1)
    if type(D).__name__=="ndarray":
        D[D>1.]=1.
        D[D<-1.]=-1.
    else:
        if D>1.: D=1.
        if D<-1.: D=-1.
    return AC(D)


class ClassVisServer():
    def __init__(self,MSName,
                 ColName="DATA",
                 TChunkSize=1,
                 TVisSizeMin=1,
                 DicoSelectOptions={},
                 LofarBeam=None,
                 AddNoiseJy=None,
                 #IdSharedMem="",
                 SM=None,NCPU=None,
                 Robust=2,Weighting="Natural",
                 WeightUVMinMax=None, WTUV=1.0,
                 GD=None,GDImag=None,
                 Name=""):

        self.GD=GD
        self.GDImag=GDImag
        self.CalcGridBasedFlags=False
        #self.IdSharedMem=IdSharedMem
        #PrefixShared="%sSharedVis"%self.IdSharedMem
        self.AddNoiseJy=AddNoiseJy
        self.ReInitChunkCount()
        self.TMemChunkSize=TChunkSize
        self.TVisSizeMin=TVisSizeMin
        self.MSName=MSName
        self.SM=SM
        self.NCPU=NCPU
        self.VisWeights=None
        self.CountPickle=0

        #global log
        #log=logger.getLogger("[%s]"%Name)
        self.Name=Name
        self.DicoSemaphores=shared_dict.create("%sDicoSemaphores"%self.Name)
        SubColName=None
        if "-" in ColName:
            ColName,SubColName=ColName.split("-")
        self.ColName=ColName
        self.SubColName=SubColName
        
        self.DicoSelectOptions=DicoSelectOptions
        # self.SharedNames=[]
        # self.PrefixShared=PrefixShared
        # self.VisInSharedMem = (PrefixShared!=None)
        self.LofarBeam=LofarBeam
        self.ApplyBeam=False
        self.DicoClusterDirs_Descriptor=None
        self.Robust=Robust
        self.Weighting=Weighting
        self.DomainsMachine=ClassJonesDomains.ClassJonesDomains()
        self.PreApplyTimesT0=np.array([],np.float64)
        self.PreApplyTimesT1=np.array([],np.float64)
        self.BeamTimesT0=np.array([],np.float64)
        self.BeamTimesT1=np.array([],np.float64)
        self.Init()
        self.dTimesVisMin=self.TVisSizeMin
        self.CurrentVisTimes_SinceStart_Sec=0.,0.

        self.iSliceSolveThis=-1


        
        self.WeightUVMinMax=WeightUVMinMax
        self.WTUV=WTUV
        # self.LoadNextVisChunk()

        # self.TEST_TLIST=[]

    def setSM(self,SM):
        self.SM=SM
        rac,decc=self.MS.radec
        if self.SM.Type=="Catalog" or self.SM.Type=="Hybrid":
            self.SM.Calc_LM(rac,decc)
            
  


            
            # if self.GD!=None:
            #     if self.GD["PreApply"]["PreApplySols"][0]!="":
            #         CJ=ClassReCluster.ClassReCluster(self.GD)
            #         CJ.ReClusterSkyModel(self.SM,self.MS.MSName)
            

    # def SetBeam(self,LofarBeam):
    #     self.BeamMode,self.DtBeamMin,self.BeamRAs,self.BeamDECs = LofarBeam
    #     useArrayFactor=("A" in self.BeamMode)
    #     useElementBeam=("E" in self.BeamMode)
    #     self.MS.LoadSR(useElementBeam=useElementBeam,useArrayFactor=useArrayFactor)
    #     self.ApplyBeam=True
    #     stop

    def Init(self,PointingID=0,NChanJones=1):
        #MSName=self.MDC.giveMS(PointingID).MSName
        kwargs={}
        DecorrMode=""
        ReadUVWDT=False
        if self.GD!=None:
            kwargs["Field"]=self.GD["DataSelection"]["FieldID"]
            kwargs["ChanSlice"]=self.GD["DataSelection"]["ChanSlice"]
            kwargs["DDID"]=self.GD["DataSelection"]["DDID"]
            DecorrMode=self.GD["SkyModel"]["Decorrelation"]
            ReadUVWDT=(("T" in DecorrMode) or ("F" in DecorrMode))
            
        self.ReadUVWDT=ReadUVWDT
        
        ToRaDec=None
        if self.GD is not None and "GDImage" in list(self.GD.keys()):
            ToRaDec=self.GD["GDImage"]["Image"]["PhaseCenterRADEC"]
        
        if ToRaDec=="align":
            log.print(ModColor.Str("kMS does not understand align mode for PhaseCenterRADEC, setting to None..."))
            ToRaDec=None
            #raise RuntimeError("incorrect BeamAt setting: use Facet or Tessel")

        if self.GD is None or not self.GD["VisData"]["ConcatFreq"]:
            from .ClassMS import ClassMS
        else:
            from .ClassMSConcat import ClassMSConcat as ClassMS
            
        MS=ClassMS(self.MSName,
                   Col=self.ColName,
                   SubCol=self.SubColName,
                   DoReadData=False,
                   ReadUVWDT=ReadUVWDT,
                   GD=self.GD,ToRADEC=ToRaDec,**kwargs)
        self.MS=MS
        self.computeSolutionTimeDomains()
        
        self.DicoMergeStations={}
        
        if self.GD and self.GD["Compression"]["MergeStations"] is not None:
            MergeStations=self.GD["Compression"]["MergeStations"]
            ListMergeNames=[]
            ListMergeStations=[]
            for Name in MergeStations:
                for iAnt in range(MS.na):
                    if Name in MS.StationNames[iAnt]:
                        ListMergeStations.append(iAnt)
                        ListMergeNames.append(MS.StationNames[iAnt])
                        
            log.print("Merging into a single station %s"%str(ListMergeNames))
            self.DicoMergeStations["ListMergeNames"]=ListMergeNames
            self.DicoMergeStations["ListBLMerge"]=[(a0,a1) for a0 in ListMergeStations for a1 in ListMergeStations if a0!=a1]
            self.DicoMergeStations["ListMergeStations"]=ListMergeStations
                        


        ######################################################
        ## Taken from ClassLOFARBeam in DDFacet
        
        ChanWidth=abs(self.MS.ChanWidth.ravel()[0])
        ChanFreqs=self.MS.ChanFreq.flatten()
        if self.GD!=None:
            NChanJones=self.GD["Solvers"]["NChanSols"]
        if NChanJones==0:
            NChanJones=self.MS.NSPWChan
        ChanEdges=np.linspace(ChanFreqs.min()-ChanWidth/2.,ChanFreqs.max()+ChanWidth/2.,NChanJones+1)

        FreqDomains=[[ChanEdges[iF],ChanEdges[iF+1]] for iF in range(NChanJones)]
        FreqDomains=np.array(FreqDomains)
        self.SolsFreqDomains=FreqDomains
        self.NChanJones=NChanJones

        MeanFreqJonesChan=(FreqDomains[:,0]+FreqDomains[:,1])/2.
        log.print("Center of frequency domains [MHz]: %s"%str((MeanFreqJonesChan/1e6).tolist()))
        DFreq=np.abs(self.MS.ChanFreq.reshape((self.MS.NSPWChan,1))-MeanFreqJonesChan.reshape((1,NChanJones)))
        self.VisToSolsChanMapping=np.argmin(DFreq,axis=1)
        log.print(("VisToSolsChanMapping %s"%DDFacet.Other.PrintList.ListToStr(self.VisToSolsChanMapping)))


        if NChanJones>self.MS.ChanFreq.size:
            raise RuntimeError("Your NChanSols [%i] is larger than the number of channels [%i] in your visibility set."%(NChanJones,self.MS.ChanFreq.size))
        
        self.SolsToVisChanMapping=[]
        for iChanSol in range(NChanJones):
            ind=np.where(self.VisToSolsChanMapping==iChanSol)[0] 
            self.SolsToVisChanMapping.append((ind[0],ind[-1]+1))
        log.print(("SolsToVisChanMapping %s"%DDFacet.Other.PrintList.ListToStr(self.SolsToVisChanMapping)))
        

        # ChanDegrid
        FreqBands=ChanEdges
        self.FreqBandsMean=(FreqBands[0:-1]+FreqBands[1::])/2.
        self.FreqBandsMin=FreqBands[0:-1].copy()
        self.FreqBandsMax=FreqBands[1::].copy()
        
        NChanDegrid = NChanJones
        MS=self.MS
        ChanDegridding=np.linspace(FreqBands.min(),FreqBands.max(),NChanDegrid+1)
        FreqChanDegridding=(ChanDegridding[1::]+ChanDegridding[0:-1])/2.
        self.FreqChanDegridding=FreqChanDegridding
        NChanDegrid=FreqChanDegridding.size
        NChanMS=MS.ChanFreq.size
        DChan=np.abs(MS.ChanFreq.reshape((NChanMS,1))-FreqChanDegridding.reshape((1,NChanDegrid)))
        ThisMappingDegrid=np.argmin(DChan,axis=1)
        self.MappingDegrid=ThisMappingDegrid
        #log.print("Mapping degrid: %s"%ListToStr(self.MappingDegrid))
        #NpShared.ToShared("%sMappingDegrid"%self.IdSharedMem,self.MappingDegrid)

        ######################################################
        
        
        # self.CalcWeigths()

        #TimesVisMin=np.arange(0,MS.DTh*60.,self.TVisSizeMin).tolist()
        #if not(MS.DTh*60. in TimesVisMin): TimesVisMin.append(MS.DTh*60.)
        #self.TimesVisMin=np.array(TimesVisMin)


    def computeSolutionTimeDomains(self):

        MS=self.MS
        print(MS)
        TStart=self.MS.F_times_all.min()
        times_0=self.MS.F_times_all-TStart
        times=np.sort(np.unique(times_0))


        TimesInt=np.arange(0,MS.DTh,self.TMemChunkSize).tolist()
        if not(MS.DTh+1./3600 in TimesInt): TimesInt.append(MS.DTh+1./3600)
        self.TimesChunk_s=np.round(np.array(TimesInt)*3600)
        TimesInt=np.array(self.TimesChunk_s)/3600
        self.TimesInt=TimesInt
        self.NTChunk=self.TimesInt.size-1
        self.TimeMemChunkRange_sec_Since70=self.TimesChunk_s+self.MS.F_times_all.min()
        L=[]
        t0=0
        TMaxObs=times.max()
        
        DicoDataDesc=OrderedDict()
        iSlice=0
        def appendSlice():
            DicoDataDesc[iSlice]={"iChunk":iChunk,
                                  "isFirstSliceOfChunk":isFirstSliceOfChunk,
                                  "t0Chunk":T0,
                                  "t1Chunk":T1,
                                  "iSlice":iSlice,
                                  "t0Slice":t0,
                                  "t1Slice":t1,
                                  "t0Chunk_70":T0+TStart,
                                  "t1Chunk_70":T1+TStart,
                                  "t0Slice_70":t0+TStart,
                                  "t1Slice_70":t1+TStart,
                                  "indRowsThisChunk":indRow}
        
        for iChunk in range(self.TimesChunk_s.size-1):
            T0=self.TimesChunk_s[iChunk]
            T1=self.TimesChunk_s[iChunk+1]
            #print("NEXT CHUNK",T0,T1)
            CChunk=((times>=T0)&(times<T1))
            ind=np.where(CChunk)[0]
            isFirstSliceOfChunk=True
            if ind.size==0:
                #print("EMPTY CHUNK")
                t0=T1
                continue
            times_s=times_0[(times_0>=T0)&(times_0<T1)]
            while True:
                t1=t0+60.*self.TVisSizeMin
                if t0>T1:
                    #print("END OF Chunk")
                    t0=T1
                    break
                if t1>T1:
                    t1=T1
                    L.append([t0,t1])
                    #print(t0,t1,ind.size,t1-t0)
                    #print("END OF Chunk")
                    #CSlice=((times_s>=t0)&(times_s<t1))
                    #ind=np.where(CSlice)[0]
                    indRow=np.where((times_s>=t0) & (times_s<t1))[0]
                    if indRow.size==0:
                        t0=T1
                        break
                    appendSlice()
                    isFirstSliceOfChunk=False
                    iSlice+=1
                    t0=T1
                    break
                indRow=np.where((times_s>=t0) & (times_s<t1))[0]
                if indRow.size==0:
                    #print("EMPTY SLICE")
                    t0=t1
                    continue
                #ts=times[ind]
                #CSlice=((times_s>=t0)&(times_s<t1))
                #ind=np.where(CSlice)[0]
                appendSlice()
                isFirstSliceOfChunk=False
                iSlice+=1
                #print(t0,t1,ind.size,t1-t0)
                t0=t1
        self.DicoDataDesc=DicoDataDesc

        iSliceLast=len(self.DicoDataDesc)-1
        D=self.DicoDataDesc[iSliceLast]
        D["t1Slice"]+=1
        D["t1Slice_70"]+=1
        D["indRowsThisChunk"]=np.where((times_s>=D["t0Slice"]) & (times_s<D["t1Slice"]))[0]
        self.NSlicePerChunk=np.zeros((self.NTChunk,),int)
        self.DicoSemaphores["SemLoadChunk"]=np.zeros((self.NTChunk,),bool)

        for iSlice in DicoDataDesc.keys():
            #print("[%.3i, %.2i] %4.1f -> %4.1f "%(iSlice,DicoDataDesc[iSlice]["iChunk"],DicoDataDesc[iSlice]["t0Slice"],DicoDataDesc[iSlice]["t1Slice"]))
            self.NSlicePerChunk[DicoDataDesc[iSlice]["iChunk"]]+=1
        
        self.NSlicesTotal=np.sum(self.NSlicePerChunk)
        self.DicoSemaphores["SemLoadSlice"]=np.zeros((self.NSlicesTotal,),bool)
        # for iChunk in range():
        #     self.DicoSemaphores["SemLoadSlice_Chunk%i"]=np.zeros((self.NTChunk,),bool)

            
        #print(self.NSlicePerChunk)

 
    def ReInitChunkCount(self):
        self.iChunkThis=-1
        self.iSliceSolveThis=-1
        
    def GiveNextVis(self):
        self.iSliceSolveThis += 1
        self.iSliceSolveNext =  self.iSliceSolveThis+1

        DataDesc_iSliceSolveThis=self.DicoDataDesc.get(self.iSliceSolveThis,None)
        iChunkThisSlice=None
        if DataDesc_iSliceSolveThis is not None:
            iChunkThisSlice=DataDesc_iSliceSolveThis["iChunk"]

        if iChunkThisSlice!=self.iChunkThis:
            self.iSliceSolveThis-=1
            return "EndChunk"
        
        isFirstSliceOfChunk=self.DicoDataDesc[self.iSliceSolveThis]["isFirstSliceOfChunk"]

        if isFirstSliceOfChunk and (iChunkThisSlice==self.iChunkThis):
            
            NameJob="%sloadNextVis:Slice=%i"%(self.Name,self.iSliceSolveThis)
            #print("run this %s"%NameJob)
            APP.APP.runJob(NameJob,
                           self._loadNextVis,io=0,
                           args=(self.iSliceSolveThis,self.iChunkThis))#,serial=True)#,
        
        DataDesc_iSliceSolveNext=self.DicoDataDesc.get(self.iSliceSolveNext,None)
        iChunkNextSlice=None
        if DataDesc_iSliceSolveNext is not None:
            iChunkNextSlice=DataDesc_iSliceSolveNext["iChunk"]
            
        if (iChunkNextSlice==self.iChunkThis):
            NameJob="%sloadNextVis:Slice=%i"%(self.Name,self.iSliceSolveNext)
            APP.APP.runJob(NameJob,
                           self._loadNextVis,
                           io=0,
                           args=(self.iSliceSolveNext,self.iChunkThis))#,serial=True)#,

        self.CurrentVisTimes_MS_Sec=self.DicoDataDesc[self.iSliceSolveThis]["t0Slice_70"],self.DicoDataDesc[self.iSliceSolveThis]["t1Slice_70"]
        
        Message=APP.APP.awaitJobResults("%sloadNextVis:Slice=%i"%(self.Name,self.iSliceSolveThis))
        #print(Message)
        if Message=="LoadSliceOK":
            DicoName="%sDicoData_%i"%(self.Name,self.iSliceSolveThis)
            D=shared_dict.attach(DicoName)
            # print("loaded %s"%DicoName)
            # print("   :",list(D.keys()))
            # print("   :",D.path,os.path.exists(D.path))
            
            # for iAnt in range(self.MS.na):
            #     ind=np.where((D["A0"]==iAnt)|(D["A1"]==iAnt))
            #     f=D["flags"][ind,:,:]
            #     print(iAnt,np.count_nonzero(f==0)/f.size)

            return D
        else:
            return Message

    def delCurrentSliceRelatedData(self):
        APP.APP.runJob("%sDelSlice_%i"%(self.Name,self.iSliceSolveThis),
                       self._delSliceRelatedData,
                       io=0,
                       args=(self.iSliceSolveThis,))#,serial=True)#,
        

    def _delSliceRelatedData(self,iSlice):
        shared_dict.delDict("%sDicoData_%i"%(self.Name,iSlice))
        shared_dict.delDict("%sDicoData.A_*.Slice_%i"%(self.Name,iSlice))
        shared_dict.delDict("%sKernelMat.A_*.Slice_%i"%(self.Name,iSlice))
        
    
    def _loadNextVis(self,iSlice,iChunk):
        SharedDicoName="%sDicoData_%i"%(self.Name,iSlice)
        D=shared_dict.attach(SharedDicoName)
        if len(D)>0:
            #print("%s already exists"%SharedDicoName)
            return "IsComputed"
        
        # print("LOAD NEXT TIME SLICE")
        # print("LOAD NEXT TIME SLICE")
        # log.print( "GiveNextVis")
        
        
        # if (t0_sec>=its_t1):
        #     return "EndChunk"
        # if not self.have_data:
        #     return "AllFlaggedThisTime"

        self.ThisDataChunk=shared_dict.attach("%sDataChunk_Chunk%i"%(self.Name,iChunk))
        T0,T1=self.ThisDataChunk["T0T1"]
        D=self.ThisDataChunk
        #if len(D)==0: stop
        
        Tmax=T1
        t0_sec=self.DicoDataDesc[iSlice]["t0Slice_70"]
        t1_sec=self.DicoDataDesc[iSlice]["t1Slice_70"]
        #print(t0_sec-Tmax)

        # time selection
        #indRowsThisChunk=self.DicoDataDesc[iSlice]["indRowsThisChunk"]
        indRowsThisChunk=np.where((self.ThisDataChunk["times"]>=t0_sec)&(self.ThisDataChunk["times"]<t1_sec))[0]

        if not np.allclose(indRowsThisChunk,self.DicoDataDesc[iSlice]["indRowsThisChunk"]): stop
        
        if indRowsThisChunk.shape[0]==0:
            if t0_sec>=Tmax:
                return "EndChunk"
            else:
                return "AllFlaggedThisTime"
            

        # np.save("indRowsThisChunk.npy",indRowsThisChunk)
        # indRowsThisChunk=np.load("indRowsThisChunk.npy")
        
            
        DATA={}
        DATA["indRowsThisChunk"]=indRowsThisChunk
        for key in D.keys():
            #print(key)
            if type(D[key])!=np.ndarray: continue
            if not(key in ['times', 'A1', 'A0', 'flags', 'uvw', 'data', 'Map_VisToJones_Time', "UVW_dt",#"IndexTimesThisChunk", 
                           "W"]):             
                DATA[key]=D[key]
            else:
                DATA[key]=D[key][indRowsThisChunk]

        #############################
        ### data selection
        #############################
        flags=DATA["flags"]
        uvw=DATA["uvw"]
        data=DATA["data"]
        A0=DATA["A0"]
        A1=DATA["A1"]
        times=DATA["times"]
        W=DATA["W"]
        Map_VisToJones_Time=DATA["Map_VisToJones_Time"]
        indRowsThisChunk=DATA["indRowsThisChunk"]
        if self.ReadUVWDT: duvw_dt=DATA["UVW_dt"]

        # IndexTimesThisChunk=DATA["IndexTimesThisChunk"]

        for Field in self.DicoSelectOptions.keys():
            if Field=="UVRangeKm":
                if self.DicoSelectOptions[Field]==None: break
                d0,d1=self.DicoSelectOptions[Field]

                d0*=1e3
                d1*=1e3
                u,v,w=uvw.T
                duv=np.sqrt(u**2+v**2)
                #ind=np.where((duv<d0)|(duv>d1))[0]
                ind=np.where((duv>d0)&(duv<d1))[0]
                
                flags=flags[ind]
                data=data[ind]
                A0=A0[ind]
                A1=A1[ind]
                uvw=uvw[ind]
                times=times[ind]

                #IndexTimesThisChunk=IndexTimesThisChunk[ind]
                W=W[ind]
                Map_VisToJones_Time=Map_VisToJones_Time[ind]
                indRowsThisChunk=indRowsThisChunk[ind]
                if self.ReadUVWDT: duvw_dt=duvw_dt[ind]

        FlagAntNumber=self.ThisDataChunk["FlagAntNumber"]
        for A in FlagAntNumber:
            ind=np.where((A0!=A)&(A1!=A))[0]
            flags=flags[ind]
            data=data[ind]
            A0=A0[ind]
            A1=A1[ind]
            uvw=uvw[ind]
            times=times[ind]
            # IndexTimesThisChunk=IndexTimesThisChunk[ind]
            W=W[ind]
            Map_VisToJones_Time=Map_VisToJones_Time[ind]
            indRowsThisChunk=indRowsThisChunk[ind]
            if self.ReadUVWDT: duvw_dt=duvw_dt[ind]
        
        if self.GD["DataSelection"]["FillFactor"]!=1.:
            Mask=np.random.rand(flags.shape[0])<self.GD["DataSelection"]["FillFactor"]
            ind=np.where(Mask)[0]
            flags=flags[ind]
            data=data[ind]
            A0=A0[ind]
            A1=A1[ind]
            uvw=uvw[ind]
            times=times[ind]
            # IndexTimesThisChunk=IndexTimesThisChunk[ind]
            W=W[ind]
            Map_VisToJones_Time=Map_VisToJones_Time[ind]
            indRowsThisChunk=indRowsThisChunk[ind]
            if self.ReadUVWDT: duvw_dt=duvw_dt[ind]
            
            

        ind=np.where(A0!=A1)[0]
        flags=flags[ind,:,:]
        data=data[ind,:,:]
        A0=A0[ind]
        A1=A1[ind]
        uvw=uvw[ind,:]
        times=times[ind]
        #IndexTimesThisChunk=IndexTimesThisChunk[ind]
        W=W[ind]
        Map_VisToJones_Time=Map_VisToJones_Time[ind]
        indRowsThisChunk=indRowsThisChunk[ind]
        if self.ReadUVWDT: duvw_dt=duvw_dt[ind]
                
        DATA["flags"]=flags
        DATA["rac_decc"]=np.array([self.MS.rac,self.MS.decc])
        DATA["uvw"]=uvw
        DATA["data"]=data
        DATA["A0"]=A0
        DATA["A1"]=A1
        DATA["times"]=times
        #DATA["IndexTimesThisChunk"]=IndexTimesThisChunk
        DATA["W"]=W
        
        DATA["Map_VisToJones_Time"]=Map_VisToJones_Time
        DATA["Map_VisToJones_Freq"]=self.ThisDataChunk["Map_VisToJones_Freq"]
        DATA["indRowsThisChunk"]=indRowsThisChunk
        DATA["iSliceChunk"]=np.array([iSlice,self.iChunkThis])
        
        if self.ReadUVWDT: DATA["UVW_dt"]=duvw_dt
                                
        #DATA["UVW_dt"]=self.MS.Give_dUVW_dt(times,A0,A1)
        
        if DATA["flags"].size==0:
            return "AllFlaggedThisTime"
        fFlagged=np.count_nonzero(DATA["flags"])/float(DATA["flags"].size)
        #print fFlagged
        if fFlagged>0.9:
            # log.print( "AllFlaggedThisTime [%f%%]"%(fFlagged*100))
            return "AllFlaggedThisTime"
        
        # if fFlagged==0.:
        #    stop
        # it0=np.min(DATA["IndexTimesThisChunk"])
        # it1=np.max(DATA["IndexTimesThisChunk"])+1
        # DATA["UVW_RefAnt"]=self.ThisDataChunk["UVW_RefAnt"][it0:it1,:,:]
        # # PM=ClassPredict(NCPU=self.NCPU,IdMemShared=self.IdSharedMem)
        # # DATA["Kp"]=PM.GiveKp(DATA,self.SM)

        # self.ClearSharedMemory()
        
        DATA=shared_dict.dict_to_shm("%sDicoData_%i"%(self.Name,iSlice), DATA)
        DATA["SolsToVisChanMapping"]=self.SolsToVisChanMapping
        DATA["A0A1"]=(DATA["A0"],DATA["A1"])
        nr,nch,npol=DATA["data"].shape
        
        #DATA["PredictArray"]=np.zeros((self.SM.NDir,nr,nch,npol),DATA["data"].dtype)
        
        #if "PreApplyJones" in D.keys():
        #    shared_dict.dict_to_shm("DicoPreApplyJones_%i"%self.iCurrentVisTime, D["PreApplyJones"])
            
        #it0=np.min(DATA["IndexTimesThisChunk"])
        #it1=np.max(DATA["IndexTimesThisChunk"])+1
        #DATA["UVW_RefAnt"]=self.ThisDataChunk["UVW_RefAnt"][it0:it1,:,:]
        
        self.DicoSemaphores["SemLoadSlice"][iSlice]=1
        #print
        #print self.MS.ROW0,self.MS.ROW1
        #t0=np.min(DATA["times"])-self.MS.F_tstart
        #t1=np.max(DATA["times"])-self.MS.F_tstart
        #self.TEST_TLIST+=sorted(list(set(DATA["times"].tolist())))
        return "LoadSliceOK"


    def setGridProps(self,Cell,nx):
        self.Cell=Cell
        self.nx=nx
        self.CalcGridBasedFlags=True



    def giveDataSizeAntenna(self):
        t=table(self.MS.MSName,ack=False)
        uvw=t.getcol("UVW")
        flags=t.getcol("FLAG")
        A0,A1=t.getcol("ANTENNA1"),t.getcol("ANTENNA2")
        t.close()
        NVisPerAnt=np.zeros(self.MS.na,np.float64)
        Field="UVRangeKm"
        self.fracNVisPerAnt=np.ones_like(NVisPerAnt)
        NVis=flags[flags==0].size
        if NVis==0:
            log.print( ModColor.Str("Hummm - All the data is flagged!!!"))
            return
        
        if self.DicoSelectOptions[Field] is not None:
            d0,d1=self.DicoSelectOptions[Field]
            
            d0*=1e3
            d1*=1e3
            u,v,w=uvw.T
            duv=np.sqrt(u**2+v**2)
            #ind=np.where((duv<d0)|(duv>d1))[0]
            ind=np.where((duv>d0)&(duv<d1))[0]
            
            flags=flags[ind]
            A0=A0[ind]
            A1=A1[ind]
            uvw=uvw[ind]

            for iAnt in range(self.MS.na):
                NVisPerAnt[iAnt]=np.where((A0==iAnt)|(A1==iAnt))[0].size

            self.fracNVisPerAnt=NVisPerAnt/np.max(NVisPerAnt)
            log.print("Fraction of data per antenna for covariance estimate: %s"%str(self.fracNVisPerAnt.tolist()))


            u,v,w=uvw.T
            d=np.sqrt(u**2+v**2)
            Compactness=np.zeros((self.MS.na,),np.float32)
            for iAnt in range(self.MS.na):
                ind=np.where((A0==iAnt)|(A1==iAnt))[0]
                if ind.size==0:
                    Compactness[iAnt]=1.e-3
                    continue
                Compactness[iAnt]=np.mean(d[ind])
            self.Compactness=Compactness/np.max(Compactness)
            log.print("Compactness: %s"%str(self.Compactness.tolist()))



            
            NVisSel=flags[flags==0].size
            log.print("Total fraction of remaining data after uv-cut: %5.2f %%"%(100*NVisSel/float(NVis)))


    # ######################################################################
            

    def LoadNextVisChunk(self):
        shared_dict.delDict("%sDicoPredictedVis_Chunk%i"%(self.Name,self.iChunkThis))
        shared_dict.delDict("%sDataChunk_Chunk%i"%(self.Name,self.iChunkThis))
        shared_dict.delDict("%sDicoPreApplyJones_Chunk%i"%(self.Name,self.iChunkThis))
        self.iChunkThis+=1
        self.iChunkNext=self.iChunkThis+1
        if self.iChunkThis==self.NTChunk:
            self.ReInitChunkCount()
            return "EndOfObservation"
        
        NameJobThis="%sloadNextChunk:Chunk=%i"%(self.Name,self.iChunkThis)
        if self.iChunkThis==0:
            APP.APP.runJob(NameJobThis,
                           self._workerLoadNextVisChunk,
                           #io=1,
                           args=(self.iChunkThis,))#,serial=True)#,

        time.sleep(0.1)
        if self.iChunkNext<self.NTChunk:
            NameJobNext="%sloadNextChunk:Chunk=%i"%(self.Name,self.iChunkNext)
            APP.APP.runJob(NameJobNext,
                           self._workerLoadNextVisChunk,
                           #io=1,
                           args=(self.iChunkNext,))#,serial=True)#,
            
        Message=APP.APP.awaitJobResults(NameJobThis)
        
        if Message=="LoadOK":
            ThisDicoPreApplyJones=shared_dict.attach("%sDicoPreApplyJones_Chunk%i"%(self.Name,self.iChunkThis))
            
            self.PreApplyTimesT0=np.concatenate([self.PreApplyTimesT0,ThisDicoPreApplyJones.get("t0",np.array([]))])
            self.PreApplyTimesT1=np.concatenate([self.PreApplyTimesT1,ThisDicoPreApplyJones.get("t1",np.array([]))])

            self.BeamTimesT0=np.concatenate([self.BeamTimesT0,ThisDicoPreApplyJones.get("BeamTimesT0",np.array([]))])
            self.BeamTimesT1=np.concatenate([self.BeamTimesT1,ThisDicoPreApplyJones.get("BeamTimesT1",np.array([]))])

            
            self.ThisDataChunk=shared_dict.attach("%sDataChunk_Chunk%i"%(self.Name,self.iChunkThis))
            self.FlagAntNumber=self.ThisDataChunk["FlagAntNumber"]
            self.TimeMemChunkRange_sec=self.ThisDataChunk["TimeMemChunkRange_sec"]
            self.TimeMemChunkRange_sec_Since70=self.ThisDataChunk["TimeMemChunkRange_sec_Since70"]
        return Message

    # def LoadNextVisChunk(self):
    #     shared_dict.delDict("DicoPredictedVis_Chunk%i"%self.iChunkThis)
    #     shared_dict.delDict("DataChunk_Chunk%i"%self.iChunkThis)
    #     shared_dict.delDict("DicoPreApplyJones_Chunk%i"%self.iChunkThis)
    #     self.iChunkThis+=1
    #     self.iChunkNext=self.iChunkThis+1
    #     if self.iChunkThis==self.NTChunk:
    #         self.ReInitChunkCount()
    #         return "EndOfObservation"
        
    #     Message=self._workerLoadNextVisChunk(self.iChunkThis)
    #     if Message=="LoadOK":
    #         self.ThisDataChunk=shared_dict.attach("DataChunk_Chunk%i"%self.iChunkThis)
    #         self.FlagAntNumber=self.ThisDataChunk["FlagAntNumber"]
    #         self.TimeMemChunkRange_sec=self.ThisDataChunk["TimeMemChunkRange_sec"]
    #         self.TimeMemChunkRange_sec_Since70=self.ThisDataChunk["TimeMemChunkRange_sec_Since70"]
    #     return Message

        
    def _workerLoadNextVisChunk(self,iChunk):
        if iChunk>0:
            while self.DicoSemaphores["SemLoadChunk"][iChunk-1]==0:
                time.sleep(0.5)
                #print("wait",iChunk)
        
        MS=self.MS

        
        # bug out when we hit the buffers
        if iChunk >= self.NTChunk:
            log.print( ModColor.Str("Reached end of observations"))
            self.ReInitChunkCount()
            self.DicoSemaphores["SemLoadChunk"][iChunk]=1
            return "EndOfObservation"

        if self.NSlicePerChunk[iChunk]==0:
            log.print( ModColor.Str("this data chunk has no slice"))
            self.DicoSemaphores["SemLoadChunk"][iChunk]=1
            return "Empty"
        
        # get current chunk boundaries
        iT0,iT1=iChunk,iChunk+1
        

        ss="Reading next data chunk in [%5.2f, %5.2f] hours  "%(self.TimesInt[iT0],self.TimesInt[iT1])
        log.print( "="*len(ss))
        log.print( ModColor.Str("Load vis chunk #%i: start"%iChunk))
        log.print(ss)
        have_data = MS.ReadData(t0=self.TimesInt[iT0],t1=self.TimesInt[iT1],ReadWeight=True)
        
        if not have_data:
            # self.CurrentVisTimes_SinceStart_Sec=self.TimesInt[iT0]*3600.,self.TimesInt[iT1]*3600.
            # self.CurrentVisTimes_MS_Sec=self.TimesInt[iT0]*3600.+self.MS.F_tstart,self.TimesInt[iT1]*3600.+self.MS.F_tstart
            log.print( ModColor.Str("this data chunk is empty"))
            self.DicoSemaphores["SemLoadChunk"][iChunk]=1
            return "Empty"

        #log.print( "    Rows= [%i, %i]"%(MS.ROW0,MS.ROW1))
        #print float(MS.ROW0)/MS.nbl,float(MS.ROW1)/MS.nbl

        ###############################
        MS=self.MS

        TimeMemChunkRange_sec=self.TimesInt[iT0]*3600.,self.TimesInt[iT1]*3600.
        TimeMemChunkRange_sec_Since70=self.TimesInt[iT0]*3600.+self.MS.F_tstart,self.TimesInt[iT1]*3600.+self.MS.F_tstart
        
        #log.print(("!!!!!!!",self.TimeMemChunkRange_sec))
        times=MS.times_all
        data=MS.data
        A0=MS.A0
        A1=MS.A1
        uvw=MS.uvw
        flags=MS.flag_all
        freqs=MS.ChanFreq.flatten()
        nbl=MS.nbl
        dfreqs=MS.dFreq
        duvw_dt=MS.uvw_dt

        # if Nchan>1:
        #     DoRevertChans=(freqs.flatten()[0]>freqs.flatten()[-1])
        # if self.DoRevertChans:
        #     print ModColor.Str("  ====================== >> Revert Channel order!")
        #     wavelength_chan=wavelength_chan[0,::-1]
        #     freqs=freqs[0,::-1]
        #     self.dFreq=np.abs(self.dFreq)
        
        
        



        #flags.fill(0)

        # f=(np.random.rand(*flags.shape)>0.5)
        # flags[f]=1
        # data[flags]=1e6

        # iAFlag=12
        # ind=np.where((A0==iAFlag)|(A1==iAFlag))[0]
        # flags[ind,:,:]=1
        
        Equidistant=IsChanEquidistant(freqs)
        if freqs.size>1:
            if Equidistant:
                log.print( "Channels are equidistant, can go fast")
            else:
                log.print( ModColor.Str("Channels are not equidistant, cannot go fast"))
        
        MS=self.MS


        if self.SM.Type=="Image":
            u,v,w=uvw.T
            wmax=self.GD["GDImage"]["CF"]["wmax"]
            wmaxkm=wmax/1000.
            log.print( "Flagging baselines with w > %f km"%(wmaxkm))
            C=299792458.
            fmax=self.MS.ChanFreq.ravel()[-1]

            ind=np.where(np.abs(w)>wmax)[0]
            flags[ind,:,:]=1

            f=ind.size/float(flags.shape[0])
            log.print( "  w-Flagged %5.1f%% of the data"%(100*f))



            # data=data[ind]
            # A0=A0[ind]
            # A1=A1[ind]
            # uvw=uvw[ind]
            # times=times[ind]
            # W=W[ind]
            # MapJones=MapJones[ind]
            # indRowsThisChunk=indRowsThisChunk[ind]


        #log.print("::::!!!!!!!!!!!!!!!")
        self.ThresholdFlag=1.#0.9
        FlagAntNumber=[]


        ########################################
        
        for A in range(MS.na):
            ind=np.where((MS.A0==A)|(MS.A1==A))[0]
            fA=MS.flag_all[ind].ravel()
            if ind.size==0:
                log.print( "Antenna #%2.2i[%s] is not in the MS"%(A,MS.StationNames[A]))
                FlagAntNumber.append(A)
                continue
                
            nf=np.count_nonzero(fA)
            
            Frac=nf/float(fA.size)
            if Frac>self.ThresholdFlag:
                log.print( "Taking antenna #%2.2i[%s] out of the solve (~%4.1f%% of flagged data, more than %4.1f%%)"%\
                    (A,MS.StationNames[A],Frac*100,self.ThresholdFlag*100))
                FlagAntNumber.append(A)
                
        if self.CalcGridBasedFlags:
            Cell=self.Cell
            nx=self.nx
            MS=self.MS
            u,v,w=MS.uvw.T
            try:
                CellRad_x,CellRad_y=(Cell/3600.)*np.pi/180
            except:
                CellRad_x=CellRad_y=(Cell/3600.)*np.pi/180
                
            d=np.sqrt((u*CellRad_x)**2+(v*CellRad_y)**2)
            #_,_,nx,ny=GridShape

            # ###
            # S=CellRad*nx
            # C=3e8
            # freqs=MS.ChanFreq
            # x=d.reshape((d.size,1))*(freqs.reshape((1,freqs.size))/C)*S
            # fA_all=(x>(nx/2))
            # ###
            
            C=3e8
            freqs=MS.ChanFreq.flatten()
            x=d.reshape((d.size,1))*(freqs.reshape((1,freqs.size))/C)
            fA_all=(x>(1./2))
            
            for A in range(MS.na):
                ind=np.where((MS.A0==A)|(MS.A1==A))[0]
                fA=fA_all[ind].ravel()
                nf=np.count_nonzero(fA)
                if fA.size==0:
                    Frac=1.0
                else:
                    Frac=nf/float(fA.size)
                if Frac>self.ThresholdFlag:
                    log.print( "Taking antenna #%2.2i[%s] out of the solve (~%4.1f%% of out-grid data, more than %4.1f%%)"%\
                               (A,MS.StationNames[A],Frac*100,self.ThresholdFlag*100))
                    FlagAntNumber.append(A)




        if self.DicoSelectOptions.get("FlagAnts",None) is not None:
            FlagAnts=self.DicoSelectOptions["FlagAnts"]
            if not isinstance(FlagAnts, list):
                FlagAnts=[FlagAnts]
            
            for Name in FlagAnts:
                for iAnt in range(MS.na):
                    if Name in MS.StationNames[iAnt]:
                        log.print( "Taking antenna #%2.2i[%s] out of the solve"%(iAnt,MS.StationNames[iAnt]))
                        FlagAntNumber.append(iAnt)



                        
        if "DistMaxToCore" in self.DicoSelectOptions.keys():
            DMax=self.DicoSelectOptions["DistMaxToCore"]*1e3
            X,Y,Z=MS.StationPos.T
            Xm,Ym,Zm=np.median(MS.StationPos,axis=0).flatten().tolist()
            D=np.sqrt((X-Xm)**2+(Y-Ym)**2+(Z-Zm)**2)
            ind=np.where(D>DMax)[0]

            for iAnt in ind.tolist():
                log.print("Taking antenna #%2.2i[%s] out of the solve (distance to core: %.1f km)"%(iAnt,MS.StationNames[iAnt],D[iAnt]/1e3))
                FlagAntNumber.append(iAnt)
            



        #############################
        #############################

        ind=np.where(np.isnan(data))
        flags[ind]=1
        
        # ## debug
        # ind=np.where((A0==0)&(A1==1))[0]
        # flags=flags[ind]
        # data=data[ind]
        # A0=A0[ind]
        # A1=A1[ind]
        # uvw=uvw[ind]
        # times=times[ind]
        # ##


        if self.AddNoiseJy!=None:
            data+=(self.AddNoiseJy/np.sqrt(2.))*(np.random.randn(*data.shape)+1j*np.random.randn(*data.shape))
            stop
        # Building uvw infos
        #################################################
        #################################################
        # log.print( "Building uvw infos .... ")
        # Luvw=np.zeros((MS.times.size,MS.na,3),uvw.dtype)
        # AntRef=0
        # indexTimes=np.zeros((times.size,),np.int64)
        # iTime=0

        # # UVW per antenna
        # #Times_all_32=np.float32(times-MS.times[0])
        # #Times32=np.float32(MS.times-MS.times[0])
        # irow=0
        # for ThisTime in MS.times:#Times32:

        #     T= ClassTimeIt.ClassTimeIt("VS")
        #     T.disable()
        #     #ind=np.where(Times_all_32[irow::]==ThisTime)[0]
        #     ind=np.where(times[irow::]==ThisTime)[0]+irow
            
        #     irow+=ind.size
        #     T.timeit("0b")

        #     indAnt=np.where(A0[ind]==AntRef)[0]
        #     ThisUVW0=uvw[ind][indAnt].copy()
        #     Ant0=A1[ind][indAnt].copy()
        #     T.timeit("1")

        #     indAnt=np.where(A1[ind]==AntRef)[0]
        #     ThisUVW1=-uvw[ind][indAnt].copy()
        #     Ant1=A0[ind][indAnt].copy()
        #     ThisUVW=np.concatenate((ThisUVW1,ThisUVW0[1::]))

        #     T.timeit("2")

        #     #AA=np.concatenate((Ant1,Ant0[1::]))
        #     Luvw[iTime,:,:]=ThisUVW[:,:]
        
        #     T.timeit("3")

        #     #Luvw.append(ThisUVW)
        #     indexTimes[ind]=iTime
        #     iTime+=1
        # log.print( "     .... Done ")
        #################################################
        #################################################


        # Dt_UVW_dt=1.*3600
        # t0=times[0]
        # t1=times[-1]
        # All_UVW_dt=[]
        # Times=np.arange(t0,t1,Dt_UVW_dt).tolist()
        # if not(t1 in Times): Times.append(t1+1)

        # All_UVW_dt=np.array([],np.float32).reshape((0,3))
        # for it in range(len(Times)-1):
        #     t0=Times[it]
        #     t1=Times[it+1]
        #     tt=(t0+t1)/2.
        #     indRows=np.where((times>=t0)&(times<t1))[0]
        #     All_UVW_dt=np.concatenate((All_UVW_dt,self.MS.Give_dUVW_dt(tt,A0[indRows],A1[indRows])))

            

        # UVW_dt=All_UVW_dt

        DoCompress=False
        if self.GD is not None:
            if (self.GD["Compression"]["CompressionMode"] is not None) or self.GD["Compression"]["CompressionDirFile"]:
                DoCompress=True

            DicoPredictedVis=None
            FreePredictColName=self.GD["VisData"]["FreePredictColName"]
            if not DoCompress and (FreePredictColName!=None):
                PredictedData=np.zeros_like(data)
                Indices=np.arange(PredictedData.size).reshape(PredictedData.shape)
                if DicoPredictedVis is None:
                    DicoPredictedVis=shared_dict.create("%sDicoPredictedVis_Chunk%i"%(self.Name,iChunk))
                DicoPredictedVis["PredictedData"]=PredictedData
                DicoPredictedVis["IndicesData"]=Indices
        
            FreePredictGainColName=self.GD["VisData"]["FreePredictGainColName"]
            if not DoCompress and (FreePredictGainColName!=None):
                PredictedDataGains=np.zeros_like(data)
                IndicesGains=np.arange(PredictedDataGains.size).reshape(PredictedDataGains.shape)
                if DicoPredictedVis is None:
                    DicoPredictedVis=shared_dict.create("%sDicoPredictedVis_Chunk%i"%(self.Name,iChunk))
                DicoPredictedVis["PredictedDataGains"]=PredictedDataGains
                DicoPredictedVis["IndicesDataGains"]=IndicesGains
        
        
        
        ThisDataChunk={"times":times,
                       "freqs":freqs,
                       "dfreqs":dfreqs,
                       "freqs_full":freqs,
                       "dfreqs_full":dfreqs,
                       #"A0A1":(A0[ind],A1[ind]),
                       #"A0A1":(A0,A1),
                       "A0":A0,
                       "A1":A1,
                       "uvw":uvw,
                       "flags":flags,
                       "nbl":nbl,
                       "na":MS.na,
                       "data":data,
                       "ROW0":MS.ROW0,
                       "ROW1":MS.ROW1,
                       "infos":np.array([MS.na,MS.TimeInterVal[0]]),
                       #"IndexTimesThisChunk":indexTimes,
                       #"UVW_RefAnt": Luvw,
                       "W":self.VisWeights[MS.ROW0:MS.ROW1],
                       #"IndRows_All_UVW_dt":IndRows_All_UVW_dt,
                       "UVW_dt":duvw_dt,
                       "T0T1":[times.min(),times.max()],
                       "FlagAntNumber":FlagAntNumber,
                       "TimeMemChunkRange_sec":TimeMemChunkRange_sec,
                       "TimeMemChunkRange_sec_Since70":TimeMemChunkRange_sec_Since70
                     }



        self.ThisDataChunk=shared_dict.dict_to_shm("%sDataChunk_Chunk%i"%(self.Name,iChunk),ThisDataChunk)
        ThisDataChunk=self.ThisDataChunk

        #self.UpdateCompression()
        #self.ThisDataChunk["Map_VisToJones_Time"]=np.zeros(([],),np.int32)
        self.ThisDataChunk["Map_VisToJones_Time"]=np.zeros((times.size,),np.int32)
        self.ThisDataChunk["Map_VisToJones_Freq"]=np.zeros((data.shape[1],),np.int32)
        self.ThisDataChunk["FlagAntNumber"]=FlagAntNumber

        ListDicoPreApply=[]
        DoPreApplyJones=False
        if self.GD is not None and self.GD["Beam"]["BeamModel"] is None:
            self.GD["Beam"]["BeamAt"]="tessel"

        #if self.SM.Type!="Image":
        #    self.GD["Beam"]["BeamAt"]="tessel"
            
            
        if self.GD is not None:

            ForceNotApply=None
            if self.SM.Type=="Catalog":
                if self.GD["Beam"]["BeamAt"].lower() == "tessel":
                    log.print("Estimating PreApply directions at the center of the tesselated areas")
                    RA,DEC=self.SM.ClusterCat.ra,self.SM.ClusterCat.dec
                    NDirJonesPreApply=RA.size
                elif self.GD["Beam"]["BeamAt"].lower() == "meansource":
                    log.print("Estimating PreApply directions at the mean sources of the tesselated areas")
                    self.SM.ComputeClusterCatWeightedPos()
                    RA,DEC=self.SM.ClusterCatMeanSource.ra,self.SM.ClusterCatMeanSource.dec
                    NDirJonesPreApply=RA.size
                elif self.GD["Beam"]["BeamAt"].lower() == "facet":
                    if self.SM.DDF_GD is None:
                        raise RuntimeError("I can't have access to the DicoImager")
                    DicoImagerFile=self.SM.DDF_GD["Output"]["Name"]
                    log.print("Estimating PreApply directions from facet infos the SM was constructed from %s.DicoFacet"%DicoImagerFile)
                    DicoImager=DDFacet.Other.MyPickle.Load("%s.DicoFacet"%DicoImagerFile)
                    RAd=np.array([DicoImager[iFacet]["RaDec"][0] for iFacet in range(len(DicoImager))])
                    DECd=np.array([DicoImager[iFacet]["RaDec"][1] for iFacet in range(len(DicoImager))])
                    
                    RA,DEC=self.SM.ClusterCat.ra,self.SM.ClusterCat.dec

                    d=AngDist(RA.reshape((-1,1)),RAd.reshape((1,-1)),
                              DEC.reshape((-1,1)),DECd.reshape((1,-1)))
                    indDirFacet=np.argmin(d,axis=1)
                    RA=RAd[indDirFacet]
                    DEC=DECd[indDirFacet]
                    NDirJonesPreApply=self.SM.ClusterCat.ra.size
                else:
                    raise RuntimeError("read the code")
            elif self.SM.Type=="Image":
                if self.GD["Beam"]["BeamAt"].lower() == "tessel":
                    log.print("Estimating PreApply directions at the center of the tesselated areas")
                    RA,DEC=self.SM.ClusterCat.ra,self.SM.ClusterCat.dec
                    NDirJonesPreApply=RA.size
                elif self.GD["Beam"]["BeamAt"].lower() == "facet":
                    log.print("Estimating PreApply directions at the center of the individual facets areas")
                    RA=np.array([self.SM.DicoImager[iFacet]["RaDec"][0] for iFacet in range(len(self.SM.DicoImager))])
                    DEC=np.array([self.SM.DicoImager[iFacet]["RaDec"][1] for iFacet in range(len(self.SM.DicoImager))])
                    NDirJonesPreApply=len(self.SM.DicoImager)
                else:
                    raise RuntimeError("read the code")
            elif self.SM.Type=="Hybrid":
                if self.GD["Beam"]["BeamAt"].lower() == "tessel":
                    log.print("Estimating PreApply directions at the center of the tesselated areas")
                    RA,DEC=self.SM.ClusterCat.ra,self.SM.ClusterCat.dec
                    NDirJonesPreApply=RA.size
                    Name=self.SM.ClusterCat.Name
                elif self.GD["Beam"]["BeamAt"].lower() == "facet":
                    log.print("Estimating PreApply directions at the center of the individual facets areas")
                    SM=self.SM.LSM[0]
                    RA=np.array([SM.DicoImager[iFacet]["RaDec"][0] for iFacet in range(len(SM.DicoImager))])
                    DEC=np.array([SM.DicoImager[iFacet]["RaDec"][1] for iFacet in range(len(SM.DicoImager))])

                    Name=SM.ClusterCat.Name
                    RA=np.concatenate([RA,self.SM.LSM[1].ClusterCat.ra])
                    DEC=np.concatenate([DEC,self.SM.LSM[1].ClusterCat.dec])
                    Name=np.concatenate([Name,self.SM.LSM[1].ClusterCat.Name])
                    NDirJonesPreApply=RA.size#len(SM.DicoImager)
                else:
                    raise RuntimeError("read the code")
                #ForceNotApply=["ATeam" in Name[i].decode("ascii") for i in range(Name.size)]
            elif self.SM.Type=="Column":
                    log.print("Estimating PreApply directions at the pointing center")
                    RA,DEC=self.MS.rac,self.MS.decc
                    NDirJonesPreApply=RA.size
                    Name=[0]
            else:
                raise RuntimeError("Incorrect option settings: SM.Type is %s, BeamAt is %s" % (self.SM.Type, self.GD["Beam"]["BeamAt"]))
            #self.SM.setDicoJonesDirToPreApplyDirs((RA,DEC))
            
            # print("[3.83185357 3.79931271 3.78953385 3.77399329] [0.57456745 0.56505037 0.56987484 0.59887151]")
            # print(RA,DEC)
            # print("3.8052541016606365 0.598298867583656")
            #################
            # if self.GD["Beam"]["BeamAt"].lower() == "tessel":
            #     log.print("Estimating PreApply directions at the center of the tesselated areas")
            #     RA,DEC=self.SM.ClusterCat.ra,self.SM.ClusterCat.dec
            #     NDirJonesPreApply=RA.size
            # # elif self.GD["Beam"]["BeamAt"].lower() == "meansource":
            # #     log.print("Estimating PreApply directions at the mean sources of the tesselated areas")
            # #     self.SM.ComputeClusterCatWeightedPos()
            # #     RA,DEC=self.SM.ClusterCatMeanSource.ra,self.SM.ClusterCatMeanSource.dec
            # #     NDirJonesPreApply=RA.size
            # elif self.GD["Beam"]["BeamAt"].lower() == "facet":
            #     if self.SM.Type!="Image":
            #         raise RuntimeError("SM must be of Image type in facet mode")
            #     log.print("Estimating PreApply directions at the center of the individual facets areas")
            #     RA=np.array([self.SM.DicoImager[iFacet]["RaDec"][0] for iFacet in range(len(self.SM.DicoImager))])
            #     DEC=np.array([self.SM.DicoImager[iFacet]["RaDec"][1] for iFacet in range(len(self.SM.DicoImager))])
            #     NDirJonesPreApply=len(self.SM.DicoImager)
            # else:
            #     raise RuntimeError("incorrect BeamAt setting: use Facet or Tessel")

            log.print("PreApply Jones matrix has %i directions"%NDirJonesPreApply)
            
            # Mapping_DirVisPredict_To_PreApplyJonesDir
            DicoBeam=None
            if self.GD["Beam"]["BeamModel"] is not None:
                
                if self.GD["Beam"]["BeamModel"]=="LOFAR":
                    NDir=RA.size
                    self.DtBeamMin=self.GD["Beam"]["DtBeamMin"]
                    useArrayFactor=("A" in self.GD["Beam"]["PhasedArrayMode"])
                    useElementBeam=("E" in self.GD["Beam"]["PhasedArrayMode"])
                    self.MS.LoadSR(useElementBeam=useElementBeam,useArrayFactor=useArrayFactor)
                    log.print( "Update LOFAR beam in %i directions [Dt = %3.1f min] ... "%(NDir,self.DtBeamMin))
                    DtBeamSec=self.DtBeamMin*60
                    tmin,tmax=np.min(times)-MS.dt/2.,np.max(times)+MS.dt/2.
                    # TimesBeam=np.arange(np.min(times),np.max(times),DtBeamSec).tolist()
                    # if not(tmax in TimesBeam): TimesBeam.append(tmax)
                    NTimesBeam=round((tmax-tmin)/DtBeamSec)
                    NTimesBeam=int(np.max([2,NTimesBeam]))


                    TimesBeam=np.linspace(np.min(times)-1,np.max(times)+1,NTimesBeam).tolist()
                    TimesBeam=np.array(TimesBeam)
                    
                    T0s=TimesBeam[:-1]
                    T1s=TimesBeam[1:]
                    Tm=(T0s+T1s)/2.
                    
                    self.BeamTimesT0=T0s
                    self.BeamTimesT1=T1s
                    # print "!!!!!!!!!!!!!!!!!!!!"
                    # T0s=MS.F_times-MS.dt/2.
                    # T1s=MS.F_times+MS.dt/2.
                    # Tm=MS.F_times

                    # from killMS.Other.rad2hmsdms import rad2hmsdms
                    # for i in range(RA.size): 
                    #     ra,dec=RA[i],DEC[i]
                    #     print rad2hmsdms(ra,Type="ra").replace(" ",":"),rad2hmsdms(dec,Type="dec").replace(" ",".")

                    Beam=np.zeros((Tm.size,NDir,self.MS.na,self.MS.NSPWChan,2,2),np.complex64)
                    for itime in range(Tm.size):
                        ThisTime=Tm[itime]
                        Beam[itime]=self.MS.GiveBeam(ThisTime,RA,DEC)
    
                    # # Beam[:,76,:,:,0,0]=20.
                    # # Beam[:,76,:,:,0,1]=0.
                    # # Beam[:,76,:,:,1,0]=0.
                    # # Beam[:,76,:,:,1,1]=20.
                    # Beam[:,:,:,:,0,0]=20.
                    # Beam[:,:,:,:,0,1]=0.
                    # Beam[:,:,:,:,1,0]=0.
                    # Beam[:,:,:,:,1,1]=20.


                    DicoBeam={}
                    DicoBeam["t0"]=T0s
                    DicoBeam["t1"]=T1s
                    DicoBeam["tm"]=Tm
                    DicoBeam["Jones"]=Beam

                    ChanWidth=self.MS.ChanWidth.ravel()[0]
                    ChanFreqs=self.MS.ChanFreq.flatten()

                    self.DomainsMachine.AddFreqDomains(DicoBeam,ChanFreqs,ChanWidth)
                    
                    NChanBeam=self.GD["Beam"]["NChanBeamPerMS"]
                    if NChanBeam==0:
                        NChanBeam=self.MS.NSPWChan
                    FreqDomainsOut=self.DomainsMachine.GiveFreqDomains(ChanFreqs,ChanWidth,NChanJones=NChanBeam)
                    self.DomainsMachine.AverageInFreq(DicoBeam,FreqDomainsOut)

                    ###### Normalise
                    rac,decc=self.MS.OriginalRadec
                    if self.GD["Beam"]["CenterNorm"]==1:
                        
                        Beam=DicoBeam["Jones"]
                        Beam0=np.zeros((Tm.size,1,self.MS.na,self.MS.NSPWChan,2,2),np.complex64)
                        for itime in range(Tm.size):
                            ThisTime=Tm[itime]
                            Beam0[itime]=self.MS.GiveBeam(ThisTime,np.array([rac]),np.array([decc]))
                            
                        DicoBeamCenter={}
                        DicoBeamCenter["t0"]=T0s
                        DicoBeamCenter["t1"]=T1s
                        DicoBeamCenter["tm"]=Tm
                        DicoBeamCenter["Jones"]=Beam0
                        self.DomainsMachine.AddFreqDomains(DicoBeamCenter,ChanFreqs,ChanWidth)
                        self.DomainsMachine.AverageInFreq(DicoBeamCenter,FreqDomainsOut)
                        Beam0=DicoBeamCenter["Jones"]
                        Beam0inv=ModLinAlg.BatchInverse(Beam0)
                        nt,nd,_,_,_,_=Beam.shape
                        Ones=np.ones((nt,nd, 1, 1, 1, 1),np.float32)
                        Beam0inv=Beam0inv*Ones
                        DicoBeam["Jones"]=ModLinAlg.BatchDot(Beam0inv,Beam)

                    ###### 




                    #nt,nd,na,nch,_,_= Beam.shape
                    #Beam=np.mean(Beam,axis=3).reshape((nt,nd,na,1,2,2))
                    
                    #DicoBeam["ChanMap"]=np.zeros((nch))
                    DicoBeam["Name"]="LOFARBeam"
                    ListDicoPreApply.append(DicoBeam)
                    
                    DoPreApplyJones=True
                    log.print( "       .... done Update LOFAR beam ")
                elif self.GD["Beam"]["BeamModel"] == "FITS" or self.GD["Beam"]["BeamModel"] == "ATCA" or self.GD["Beam"]["BeamModel"] == "NENUFAR":
                    NDir = RA.size
                    self.DtBeamMin = self.GD["Beam"]["DtBeamMin"]

                    if self.GD["Beam"]["BeamModel"] == "FITS":
                        from DDFacet.Data.ClassFITSBeam import ClassFITSBeam as ClassDDFBeam
                    elif self.GD["Beam"]["BeamModel"] == "ATCA":
                        from DDFacet.Data.ClassATCABeam import ClassATCABeam as ClassDDFBeam
                    elif self.GD["Beam"]["BeamModel"] == "NENUFAR":
                        from DDFacet.Data.ClassNenuBeam import ClassNenuBeam as ClassDDFBeam
                        
                    # make fake opts dict (DDFacet clss expects slightly different option names)
                    opts = self.GD["Beam"]
                    opts["NBand"] = opts["NChanBeamPerMS"]
                    ddfbeam = ClassDDFBeam(self.MS, opts)

                    TimesBeam = np.array(ddfbeam.getBeamSampleTimes(times))
                    FreqDomains = ddfbeam.getFreqDomains()
                    nfreq_dom = FreqDomains.shape[0]

                    log.print( "Update %s beam in %i dirs, %i times, %i freqs ... " % (self.GD["Beam"]["BeamModel"],NDir, len(TimesBeam), nfreq_dom))

                    T0s = TimesBeam[:-1]
                    T1s = TimesBeam[1:]
                    Tm = (T0s + T1s) / 2.

                    #self.BeamTimes = TimesBeam
                    self.BeamTimesT0=T0s
                    self.BeamTimesT1=T1s

                    Beam = np.zeros((Tm.size, NDir, self.MS.na, FreqDomains.shape[0], 2, 2), np.complex64)
                    for itime, tm in enumerate(Tm):
                        Beam[itime] = ddfbeam.evaluateBeam(tm, RA, DEC)

                    DicoBeam = {}
                    DicoBeam["t0"] = T0s
                    DicoBeam["t1"] = T1s
                    DicoBeam["tm"] = Tm
                    DicoBeam["Jones"] = Beam
                    DicoBeam["FreqDomain"] = FreqDomains

                    ###### Normalise
                    #rac, decc = self.MS.radec
                    rac,decc=self.MS.OriginalRadec
                    if self.GD["Beam"]["CenterNorm"] == 1:

                        Beam = DicoBeam["Jones"]
                        Beam0 = np.zeros((Tm.size, 1, self.MS.na, nfreq_dom, 2, 2), np.complex64)
                        for itime, tm in enumerate(Tm):
                            Beam0[itime] = ddfbeam.evaluateBeam(tm, np.array([rac]), np.array([decc]))

                        DicoBeamCenter = {}
                        DicoBeamCenter["t0"] = T0s
                        DicoBeamCenter["t1"] = T1s
                        DicoBeamCenter["tm"] = Tm
                        DicoBeamCenter["Jones"] = Beam0
                        DicoBeamCenter["FreqDomain"] = FreqDomains
                        Beam0inv = ModLinAlg.BatchInverse(Beam0)
                        nt, nd, _, _, _, _ = Beam.shape
                        Ones = np.ones((nt, nd, 1, 1, 1, 1), np.float32)
                        Beam0inv = Beam0inv * Ones
                        DicoBeam["Jones"] = ModLinAlg.BatchDot(Beam0inv, Beam)

                    ######




                    # nt,nd,na,nch,_,_= Beam.shape
                    # Beam=np.mean(Beam,axis=3).reshape((nt,nd,na,1,2,2))

                    # DicoBeam["ChanMap"]=np.zeros((nch))
                    ListDicoPreApply.append(DicoBeam)

                    DoPreApplyJones = True
                    log.print( "       .... done Update beam ")
                elif self.GD["Beam"]["BeamModel"] == "GMRT":
                    NDir = RA.size
                    self.DtBeamMin = self.GD["Beam"]["DtBeamMin"]

                    from DDFacet.Data.ClassGMRTBeam import ClassGMRTBeam
                    # make fake opts dict (DDFacet clss expects slightly different option names)
                    opts = self.GD["Beam"]
                    opts["NBand"] = opts["NChanBeamPerMS"]
                    gmrtbeam = ClassGMRTBeam(self.MS, opts)

                    TimesBeam = np.array(gmrtbeam.getBeamSampleTimes(times))
                    FreqDomains = gmrtbeam.getFreqDomains()
                    nfreq_dom = FreqDomains.shape[0]

                    log.print( "Update GMRT beam in %i dirs, %i times, %i freqs ... " % (NDir, len(TimesBeam), nfreq_dom))

                    T0s = TimesBeam[:-1]
                    T1s = TimesBeam[1:]
                    Tm = (T0s + T1s) / 2.

                    #self.BeamTimes = TimesBeam
                    self.BeamTimesT0=T0s
                    self.BeamTimesT1=T1s

                    Beam = np.zeros((Tm.size, NDir, self.MS.na, FreqDomains.shape[0], 2, 2), np.complex64)
                    for itime, tm in enumerate(Tm):
                        Beam[itime] = gmrtbeam.GiveInstrumentBeam(tm, RA, DEC)

                    DicoBeam = {}
                    DicoBeam["t0"] = T0s
                    DicoBeam["t1"] = T1s
                    DicoBeam["tm"] = Tm
                    DicoBeam["Jones"] = Beam
                    DicoBeam["FreqDomain"] = FreqDomains

                    ###### Normalise
                    #rac, decc = self.MS.radec
                    rac,decc=self.MS.OriginalRadec
                    if self.GD["Beam"]["CenterNorm"] == 1:

                        Beam = DicoBeam["Jones"]
                        Beam0 = np.zeros((Tm.size, 1, self.MS.na, nfreq_dom, 2, 2), np.complex64)
                        for itime, tm in enumerate(Tm):
                            Beam0[itime] = gmrtbeam.evaluateBeam(tm, np.array([rac]), np.array([decc]))

                        DicoBeamCenter = {}
                        DicoBeamCenter["t0"] = T0s
                        DicoBeamCenter["t1"] = T1s
                        DicoBeamCenter["tm"] = Tm
                        DicoBeamCenter["Jones"] = Beam0
                        DicoBeamCenter["FreqDomain"] = FreqDomains
                        Beam0inv = ModLinAlg.BatchInverse(Beam0)
                        nt, nd, _, _, _, _ = Beam.shape
                        Ones = np.ones((nt, nd, 1, 1, 1, 1), np.float32)
                        Beam0inv = Beam0inv * Ones
                        DicoBeam["Jones"] = ModLinAlg.BatchDot(Beam0inv, Beam)

                    ######




                    # nt,nd,na,nch,_,_= Beam.shape
                    # Beam=np.mean(Beam,axis=3).reshape((nt,nd,na,1,2,2))

                    # DicoBeam["ChanMap"]=np.zeros((nch))
                    ListDicoPreApply.append(DicoBeam)

                    DoPreApplyJones = True
                    log.print( "       .... done Update GMRT beam ")

                if ForceNotApply is not None:
                    for iDir,DoForceUnity in enumerate(ForceNotApply):
                        if DoForceUnity:
                            log.print("  [%s] Setting Beam Jones matrices to unity"%(Name[iDir].decode("ascii")))
                            DicoBeam["Jones"][:,iDir,:,:,:,:].fill(0)
                            DicoBeam["Jones"][:,iDir,:,:,0,0]=1
                            DicoBeam["Jones"][:,iDir,:,:,1,1]=1

                
            if isinstance(self.GD["PreApply"]["PreApplySols"],str):
                self.GD["PreApply"]["PreApplySols"]=[self.GD["PreApply"]["PreApplySols"]]
            if isinstance(self.GD["PreApply"]["PreApplyMode"],str):
                self.GD["PreApply"]["PreApplyMode"]=[self.GD["PreApply"]["PreApplyMode"]]
                
            if (self.GD["PreApply"]["PreApplySols"] is not None) and (self.GD["PreApply"]["PreApplySols"][0]!=""):
                ModeList=self.GD["PreApply"]["PreApplyMode"]
                if ModeList==[""]:
                    ModeList=["AP"]*len(self.GD["PreApply"]["PreApplySols"])
                
                for SolFile,Mode in zip(self.GD["PreApply"]["PreApplySols"],ModeList):
                    
                    # if (SolFile!="")&(not(".npz" in SolFile)):
                    #     Method=SolFile
                    #     ThisMSName=reformat.reformat(os.path.abspath(self.MS.MSName),LastSlash=False)
                    #     SolFileLoad="%s/killMS.%s.sols.npz"%(ThisMSName,Method)
                    #     if self.GD["Solutions"]["SolsDir"]:
                    #         _MSName=reformat.reformat(self.MSName).split("/")[-2]
                    #         DirName="%s%s"%(reformat.reformat(self.GD["Solutions"]["SolsDir"]),_MSName)
                    #         SolFileLoad="%s/killMS.%s.sols.npz"%(DirName,SolFile)
                    # else:
                    #     SolFileLoad=SolFile

                    CGiveSaveFileName=ClassGiveSolsFile.ClassGive_kMSFileName(GD=self.GD)
                    SolFileLoad=CGiveSaveFileName.GiveFileName(SolFile,Type="Sols")
                    log.print( "Loading solution file %s in %s mode"%(SolFileLoad,Mode))

                    S=np.load(SolFileLoad)
                    Sols=S["Sols"]
                    ThisClusterCat=S["ClusterCat"].view(np.recarray)
                    
                    if RA is not None:
                        if ThisClusterCat.ra.size==RA.size:
                            C0=np.allclose(ThisClusterCat.ra,RA)
                            C1=np.allclose(ThisClusterCat.dec,DEC)
                            RAd=ThisClusterCat.ra
                            DECd=ThisClusterCat.dec
                            d=AngDist(RA.reshape((-1,1)),RAd.reshape((1,-1)),
                                      DEC.reshape((-1,1)),DECd.reshape((1,-1)))
                            d*=np.pi/180
                            if (not C0 or not C1) and (ThisClusterCat.ra.size!=1):
                                log.print(ModColor.Str("incompatible directions for RA/DEC Jones merging (max d=%.4f deg)"%d.max()))
                            
                                
                        
                    
                    


                    
                    nt,nch,na,nd,_,_=Sols["G"].shape
                    
                    DicoSols={}
                    DicoSols["t0"]=Sols["t0"]
                    DicoSols["t1"]=Sols["t1"]
                    DicoSols["tm"]=(Sols["t0"]+Sols["t1"])/2.
                    DicoSols["Jones"]=np.swapaxes(Sols["G"],1,3).reshape((nt,nd,na,nch,2,2))

                    #NDirJonesPreApplyNDirPreApply=self.SM.NDirsOrig
                    if NDirJonesPreApply!=nd:
                        if nd==1:
                            log.print(("  PreApplySol %s is DI, adapting shape to the %i directions of the PreApply Jones"%
                                       (SolFile,NDirJonesPreApply)))
                            DicoSols["Jones"]=DicoSols["Jones"]*np.ones((1,NDirJonesPreApply,1,1,1,1))
                        else:
                            log.print("   [%i]  JSolved  = J0_preAp [%i] ([]directions)"%(self.SM.NDir,NDirJonesPreApply))
                            log.print("   [%i] J0_preAp *= J1_file  [%i] "%(NDirJonesPreApply,nd))
                            log.print("      finding the nearest directions:")
                            JJ=np.zeros((nt,NDirJonesPreApply,na,nch,2,2),DicoSols["Jones"].dtype)
                            RAd=ThisClusterCat.ra
                            DECd=ThisClusterCat.dec
                            d=AngDist(RA.reshape((-1,1)),RAd.reshape((1,-1)),
                                      DEC.reshape((-1,1)),DECd.reshape((1,-1)))
                            for iDir in range(NDirJonesPreApply):
                                iDirJones=np.argmin(d[iDir])
                                dmin=d[iDir][iDirJones]*180/np.pi
                                log.print("        J0_preAp[%i] *= J1_file[%i] (d=%.2f deg away)"%(iDir,iDirJones,dmin))
                                JJ[:,iDir,...]=DicoSols["Jones"][:,iDirJones,...]
                            DicoSols["Jones"]=JJ
                            
                            
                            
                    DicoSols["FreqDomain"]=S["FreqDomains"]
                    if not("A" in Mode):
                        log.print("  normalising Jones matrices by their amplitude")
                        ind=(DicoSols["Jones"]!=0.)
                        DicoSols["Jones"][ind]/=np.abs(DicoSols["Jones"][ind])
                        
                    if not("P" in Mode):
                        dtype=DicoSols["Jones"].dtype
                        log.print("  zeroing the phases of the Jones matrices")
                        DicoSols["Jones"]=(np.abs(DicoSols["Jones"]).astype(dtype)).copy()

                    #DicoSols["Jones"]=Sols["G"].reshape((nt,nd,na,1,2,2))
                    DicoSols["Name"]=SolFileLoad
                    ListDicoPreApply.append(DicoSols)
                    DoPreApplyJones=True


                    
        if DoPreApplyJones:
            DicoJones=ListDicoPreApply[0]
            DomainsMachine=self.DomainsMachine

            if self.SM.Type=="Image":
                DomainsMachine.setFacetToDirMapping([self.SM.DicoImager[iFacet]["iDirJones"] for iFacet in range(len(self.SM.DicoImager))])

            

            
            for DicoJones1 in ListDicoPreApply[1::]:
                DicoJones=DomainsMachine.MergeJones(DicoJones1,DicoJones)
                
            
            Map_VisToJones_Time,Map_VisToJones_Freq=DomainsMachine.giveVisToJonesMapping(DicoJones,self.ThisDataChunk["times"],self.ThisDataChunk["freqs"])
            self.ThisDataChunk["Map_VisToJones_Time"]=Map_VisToJones_Time
            self.ThisDataChunk["Map_VisToJones_Freq"]=Map_VisToJones_Freq
            
            
            log.print(("VisToSolsChanMapping %s"%DDFacet.Other.PrintList.ListToStr(Map_VisToJones_Freq)))
            
            # ind=np.zeros((times.size,),np.int32)
            # #nt,na,nd,_,_,_=Beam.shape
            # ii=0
            # for it in range(nt):
            #     t0=DicoJones["t0"][it]
            #     t1=DicoJones["t1"][it]
            #     indMStime=np.where((times>=t0)&(times<t1))[0]
            #     indMStime=np.ones((indMStime.size,),np.int32)*it
            #     ind[ii:ii+indMStime.size]=indMStime[:]
            #     ii+=indMStime.size
            # TimeMapping=ind

            #DicoJones["ChanMap"]=self.VisToJonesChanMapping
            #self.ThisDataChunk["MapJones"]=TimeMapping




            DicoClusterDirs={}
            # DicoClusterDirs["l"]=self.SM.ClusterCatOrig.l
            # DicoClusterDirs["m"]=self.SM.ClusterCatOrig.m
            # DicoClusterDirs["ra"]=self.SM.ClusterCatOrig.ra
            # DicoClusterDirs["dec"]=self.SM.ClusterCatOrig.dec
            # DicoClusterDirs["I"]=self.SM.ClusterCatOrig.SumI
            # DicoClusterDirs["Cluster"]=self.SM.ClusterCatOrig.Cluster
            DicoClusterDirs["ra"]=RA
            DicoClusterDirs["dec"]=DEC
            
            #DicoClusterDirs=shared_dict.dict_to_shm("DicoClusterDirs",DicoClusterDirs)
            
            DicoPreApplyJones=DicoJones
            DicoPreApplyJones["DicoFacetDirs"]=DicoClusterDirs

                        

            
            #self.ThisDataChunk["Map_VisToJones_Time"]=DicoJones["Map_VisToJones_Time"].copy()
            #del(DicoJones["Map_VisToJones_Time"])
            DicoPreApplyJones["BeamTimesT0"]=self.BeamTimesT0
            DicoPreApplyJones["BeamTimesT1"]=self.BeamTimesT1

            shared_dict.dict_to_shm("%sDicoPreApplyJones_Chunk%i"%(self.Name,iChunk),DicoPreApplyJones)

        log.print( ModColor.Str("Load vis chunk #%i: done"%iChunk))
            

        self.DicoSemaphores["SemLoadChunk"][iChunk]=1
        return "LoadOK"


    def setFOV(self,FullImShape,PaddedFacetShape,FacetShape,CellSizeRad):
        self.FullImShape=FullImShape
        self.PaddedFacetShape=PaddedFacetShape
        self.FacetShape=FacetShape
        self.CellSizeRad=CellSizeRad

    # ###########################

    def GiveAllUVW(self):
        if self.GD is not None:
            if self.GD["Weighting"]["WeightInCol"] is not None and self.GD["Weighting"]["WeightInCol"]!="":
                WeightCol=self.GD["Weighting"]["WeightInCol"]
            else:
                WeightCol="WEIGHT"
        else:
            WeightCol="WEIGHT"
        return self.MS.GiveAllUVW(WeightCol=WeightCol)

    def CalcWeigths(self,FOV=5.):

        if self.VisWeights!=None: return
        
        uvw,WEIGHT,flags=self.GiveAllUVW()
        u,v,w=uvw.T

        freq=np.mean(self.MS.ChanFreq)
        uvdist=np.sqrt(u**2+v**2)
        uvmax=np.max(uvdist)
        CellSizeRad=res=1./(uvmax*freq/3.e8)
        npix=(FOV*np.pi/180)/res
        
        npix=np.min([npix,30000])

        ImShape=(1,1,npix,npix)
        #VisWeights=WEIGHT[:,0]#np.ones((uvw.shape[0],),dtype=np.float32)
        VisWeights=np.ones((uvw.shape[0],),dtype=np.float32)
        Robust=self.Robust

        # WeightMachine=ClassWeighting.ClassWeighting(ImShape,res)
        # self.VisWeights=WeightMachine.CalcWeights(uvw,VisWeights,Robust=Robust,
        #                                           Weighting=self.Weighting)

        ######################
        #uvw,WEIGHT,flags=self.GiveAllUVW()

        if self.GD is not None:
            if self.GD["Weighting"]["WeightInCol"] is not None and self.GD["Weighting"]["WeightInCol"]!="":
                log.print("Using column %s to compute the weights"%self.GD["Weighting"]["WeightInCol"])
                VisWeights=WEIGHT
            else:
                VisWeights=np.ones((uvw.shape[0],self.MS.ChanFreq.size),dtype=np.float32)
        else:
            VisWeights=np.ones((uvw.shape[0],self.MS.ChanFreq.size),dtype=np.float32)
        
        if np.max(VisWeights)==0.:
            log.print("All imaging weights are 0, setting them to ones")
            VisWeights.fill(1)

        if self.SM.Type=="Image":
            ImShape=self.PaddedFacetShape
            CellSizeRad=self.CellSizeRad

        WeightMachine=ClassWeighting.ClassWeighting(ImShape,CellSizeRad)#res)
        VisWeights=WeightMachine.CalcWeights(uvw,VisWeights,flags,self.MS.ChanFreq,
                                             Robust=Robust,
                                             Weighting=self.Weighting)

        if self.WeightUVMinMax is not None:
            uvmin,uvmax=self.WeightUVMinMax
            log.print('Giving full weight to data in range %f - %f km' % (uvmin, uvmax))
            uvmin*=1000
            uvmax*=1000
            filter=(uvdist<uvmin) | (uvdist>uvmax)
            log.print('Downweighting %i out of %i visibilities' % (np.sum(filter),len(uvdist)))
            VisWeights[filter]*=self.WTUV

        MeanW=np.mean(VisWeights[VisWeights!=0.])
        VisWeights/=MeanW
        log.print('Min weight is %f max is %f' % (np.min(VisWeights),np.max(VisWeights)))
        #VisWeight[VisWeight==0.]=1.
        self.VisWeights=VisWeights
    
