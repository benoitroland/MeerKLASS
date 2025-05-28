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
from killMS.Data.ClassMS import ClassMS
#from Sky.ClassSM import ClassSM
from killMS.Other.ClassTimeIt import ClassTimeIt
import numexpr as ne
#import ModNumExpr
import scipy.special
from DDFacet.Other import AsyncProcessPool as APP

import multiprocessing
from killMS.Array import ModLinAlg

#ne.evaluate=lambda sin: ("return %s"%sin)
import time
from DDFacet.Array import shared_dict

import six
if six.PY2:
    try:
        from killMS.Predict import predict27 as predict
    except ImportError:
        from killMS.cbuild.Predict import predict27 as predict
elif six.PY3:
    from killMS.cbuild.Predict import predict3x as predict 

from killMS.Other import findrms
from killMS.Other import ModColor
from killMS.Other.ModChanEquidistant import IsChanEquidistant

try:
    from DDFacet.Array import shared_dict
    from DDFacet.Imager import ClassDDEGridMachine
except:
    pass

# def SolsToDicoJones(Sols,nf):
#     Jones={}
#     Jones["t0"]=Sols.t0
#     Jones["t1"]=Sols.t1
#     nt,na,nd,_,_=Sols.G.shape
#     G=np.swapaxes(Sols.G,1,2).reshape((nt,nd,na,1,2,2))
#     Jones["Beam"]=G
#     Jones["BeamH"]=ModLinAlg.BatchH(G)
#     Jones["ChanMap"]=np.zeros((nf,))#.tolist()
#     return Jones





####################################################
####################################################

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


class ClassPredict():
    def __init__(self,Precision="S",NCPU=6,IdMemShared=None,DoSmearing="",BeamAtFacet=False,LExp=None,LSinc=None,LBessel=None):
        self.NCPU=NCPU
        ne.set_num_threads(self.NCPU)
        if Precision=="D":
            self.CType=np.complex128
            self.FType=np.float64
        if Precision=="S":
            self.CType=np.complex64
            self.FType=np.float32
        self.DoSmearing=DoSmearing
        self.IdSharedMem=IdMemShared
        self._BeamAtFacet = BeamAtFacet
        self.SM=None
        
        Np=2
        if self.DoSmearing!=0:
            if (("F" in DoSmearing)or("T" in DoSmearing)): Np=100000

        if LBessel==None:
            x=np.linspace(1e-6,100,10000)
            y=2.*scipy.special.jv(1,x)/x
            LBessel=[np.float32(y),np.float64(x[1]-x[0])]

        self.LBessel=LBessel
            
        if LExp==None:
            x=np.linspace(0.,10,Np)
            Exp=np.float32(np.exp(-x))
            LExp=[Exp,x[1]-x[0]]
        self.LExp=LExp

        if LSinc==None:
            x=np.linspace(0.,10,Np)
            Sinc=np.zeros(x.shape,np.float32)
            Sinc[0]=1.
            Sinc[1::]=np.sin(x[1::])/(x[1::])
            LSinc=[Sinc,x[1]-x[0]]
            
            #phi=1.471034
            #d=int(phi/(x[1]-x[0]))
            #print Sinc[d],np.sin(phi)/(phi)
            #stop
        self.LSinc=LSinc

    def setSM(self,SM):
        self.SM=SM


    def ApplyCal(self,DicoData,ApplyTimeJones,iCluster):
        D=ApplyTimeJones
        Jones=D["Jones"]
        JonesH=D["JonesH"]
        lt0,lt1=D["t0"],D["t1"]
        ColOutDir=DicoData["data"]
        A0=DicoData["A0"]
        A1=DicoData["A1"]
        times=DicoData["times"]
        na=int(DicoData["infos"][0])

        
        for it in range(lt0.size):
            
            t0,t1=lt0[it],lt1[it]
            ind=np.where((times>=t0)&(times<t1))[0]

            if ind.size==0: continue
            data=ColOutDir[ind]
            # flags=DicoData["flags"][ind]
            A0sel=A0[ind]
            A1sel=A1[ind]
            
            if "Map_VisToJones_Freq" in DicoData.keys():
                ChanMap=DicoData["Map_VisToJones_Freq"]
            else:
                ChanMap=range(nf)

            for ichan in range(len(ChanMap)):
                JChan=ChanMap[ichan]
                if iCluster!=-1:
                    J0=Jones[it,iCluster,:,JChan,:,:].reshape((na,4))
                    JH0=JonesH[it,iCluster,:,JChan,:,:].reshape((na,4))
                else:
                    J0=np.mean(Jones[it,:,:,JChan,:,:],axis=1).reshape((na,4))
                    JH0=np.mean(JonesH[it,:,:,JChan,:,:],axis=1).reshape((na,4))
                
                J=ModLinAlg.BatchInverse(J0)
                JH=ModLinAlg.BatchInverse(JH0)

                data[:,ichan,:]=ModLinAlg.BatchDot(J[A0sel,:],data[:,ichan,:])
                data[:,ichan,:]=ModLinAlg.BatchDot(data[:,ichan,:],JH[A1sel,:])

                # Abs_g0=(np.abs(J0[A0sel,0])<Threshold)
                # Abs_g1=(np.abs(JH0[A1sel,0])<Threshold)
                # flags[Abs_g0,ichan,:]=True
                # flags[Abs_g1,ichan,:]=True

            ColOutDir[ind]=data[:]

            # DicoData["flags"][ind]=flags[:]


    def GiveParamJonesList(self,DicoJonesMatricesIn,A0,A1,Map_VisToJones_Time,Map_VisToJones_Freq):
        if "DicoApplyJones" in DicoJonesMatricesIn.keys():
            DicoJonesMatrices=DicoJonesMatricesIn["DicoApplyJones"]
        else:
            DicoJonesMatrices=DicoJonesMatricesIn
            
        JonesMatrices=np.complex64(DicoJonesMatrices["Jones"])
        Map_VisToJones_Time=np.int32(Map_VisToJones_Time)#DicoJonesMatrices["Map_VisToJones_Time"])
        # print DicoJonesMatrices.keys()
        # print DicoJonesMatricesIn.keys()

        Map_VisToJones_Freq=np.int32(Map_VisToJones_Freq)#DicoJonesMatrices["Map_VisToJones_Freq"])
        #MapJones=np.int32(np.arange(A0.shape[0]))
        #print DicoJonesMatrices["MapJones"].shape
        #stop
        A0=np.int32(A0)
        A1=np.int32(A1)
        if not np.allclose(A0.size,A1.size,Map_VisToJones_Time.size):
            stop
        #print Map_VisToJones_Time
        #print Map_VisToJones_Freq
        ParamJonesList=[Map_VisToJones_Time,A0,A1,JonesMatrices,Map_VisToJones_Freq]
        # print sorted(list(set(ParamJonesList[0].tolist())))
        return ParamJonesList

    def GiveCovariance(self,DicoData,ApplyTimeJones,SM):
        if self.SM is not None:
            SM=self.SM
        
        D=ApplyTimeJones
        Beam=D["Beam"]
        BeamH=D["BeamH"]
        lt0,lt1=D["t0"],D["t1"]
        A0=DicoData["A0"]
        A1=DicoData["A1"]
        times=DicoData["times"]
        na=DicoData["infos"][0]
        
        
        #print "predict...",times[0],times[-1]
        #Predict=self.predictKernelPolCluster(DicoData,SM,ApplyTimeJones=ApplyTimeJones)
        #print "done..."


        Resid=DicoData["resid"]#-Predict
        flags=DicoData["flags"]

        ListDirection=SM.Dirs

        #import pylab
        #pylab.clf()
        #import ppgplot

        MaxMat=np.zeros(Resid.shape,dtype=np.float32)
        CVis=np.zeros_like(Resid)

        


        for iCluster in ListDirection:
            CVis.fill(0)
            ParamJonesList=self.GiveParamJonesList(ApplyTimeJones,A0,A1,Map_VisToJones_Time,Map_VisToJones_Freq)
            ParamJonesList=ParamJonesList+[iCluster]
            predict.CorrVis(Resid,CVis,ParamJonesList)

            CVis[flags==1]=0.
            aCVis=np.abs(CVis)
            
            #print "In direction %i: (std, std_abs, med_abs, max_abs)=(%f, %f, %f, %f)"%(iCluster,np.std(CVis),np.std(aCVis),np.median(aCVis),np.max(aCVis))
            ind=(aCVis>MaxMat)
            MaxMat[ind]=aCVis[ind]

            del(aCVis)


        nrow,nch,_=Resid.shape
        W=DicoData["W"]#np.ones((nrow,nch),np.float32)
        rms=findrms.findrms(MaxMat)
        diff=np.abs(MaxMat-np.median(MaxMat))/rms
        cond=(diff>3.)
        ind=np.any(cond,axis=2)
        W[ind]=0.


    def GiveResidAntCovariance(self,DicoData,ApplyTimeJones,SM):
        if self.SM is not None:
            SM=self.SM
            
        D=ApplyTimeJones
        Jones=D["Jones"]
        JonesH=D["JonesH"]
        lt0,lt1=D["t0"],D["t1"]
        ColOutDir=DicoData["data"]
        A0=DicoData["A0"]
        A1=DicoData["A1"]
        times=DicoData["times"]
        na=DicoData["infos"][0]
        Sigma=ApplyTimeJones["Stats"]
        W=DicoData["W"]#np.ones((nrow,nch),np.float32)
        nrow,nch=W.shape

        # print Jones.shape
        rmsAllAnts=Sigma[:,:,:,0]


        rmsAllAnts=rmsAllAnts[rmsAllAnts>0.]
        
        if len(rmsAllAnts)==0:
            W.fill(1)
            return
        rms=np.min(rmsAllAnts)
        #print rmsAllAnts,rms
        S=Sigma[:,:,:,0]
        S[S==0]=1e6
        S[S==-1]=1e6


        for it in range(lt0.size):
            
            t0,t1=lt0[it],lt1[it]
            ind=np.where((times>=t0)&(times<t1))[0]

            if ind.size==0: continue
            data=ColOutDir[ind]
            # flags=DicoData["flags"][ind]

            # print ind,t0,t1

            A0sel=A0[ind]
            A1sel=A1[ind]
            
            ThisStats=np.sqrt(np.abs(Sigma[it][0,:,0]**2-rms**2))
            # try:
            #     ThisStats=np.sqrt(np.abs(Sigma[it][:,0]**2-rms**2))
            # except:
            #     print "S",Sigma[it][:,0]
            #     print "rms",rms
            #     print "Resid",(Sigma[it][:,0]**2-rms**2)
            #     #ApplyTimeJones["Stats"]

            # print "t:",it
            # print "0",np.min(Sigma[it][:,0]**2-rms**2)
            # print "1",np.min(Sigma[it][:,0]**2)
            # print "2",(rms**2)

            # nt,nd,na,1,2,2
            Jabs=np.abs(Jones[it,:,:,0,0,0])
            J=np.mean(Jabs,axis=0)
            # print "================="
            # print A0sel
            # print Jones.shape,Jabs.shape,J.shape
            J0=J[A0sel]
            J1=J[A1sel]
            
            sig0=ThisStats[A0sel]
            sig1=ThisStats[A1sel]

            fact=1.
            w=1./(rms**2+(sig0)**(2*fact)+(sig1)**(2*fact))#+sig0*sig1)
            for ich in range(nch):
                W[ind,ich]=w[:]




    # def GiveGM(self,iFacet,SM):
    #     GridMachine=ClassDDEGridMachine.ClassDDEGridMachine(SM.GD,
    #                                                         SM.DicoImager[iFacet]["DicoConfigGM"]["ChanFreq"],
    #                                                         SM.DicoImager[iFacet]["DicoConfigGM"]["Npix"],
    #                                                         lmShift=SM.DicoImager[iFacet]["lmShift"],
    #                                                         IdSharedMem=self.IdSharedMem,
    #                                                         IDFacet=SM.DicoImager[iFacet]["IDFacet"],
    #                                                         SpheNorm=False)

    def GiveGM(self,iFacet,SM):
        """
        Factory: Initializes a gridding machine for this facet
        Args:
            iFacet: index of facet

        Returns:
            grid machine instance
        """

        # GridMachine=ClassDDEGridMachine.ClassDDEGridMachine(SM.GD,#RaDec=self.DicoImager[iFacet]["RaDec"],
        #                                                     SM.DicoImager[iFacet]["DicoConfigGM"]["ChanFreq"],
        #                                                     SM.DicoImager[iFacet]["DicoConfigGM"]["NPix"],
        #                                                     lmShift=SM.DicoImager[iFacet]["lmShift"],
        #                                                     IdSharedMem=IdSharedMem,
        #                                                     IdSharedMemData=IdSharedMemData,
        #                                                     FacetDataCache=FacetDataCache,
        #                                                     ChunkDataCache=ChunkDataCache,
        #                                                     IDFacet=SM.DicoImager[iFacet]["IDFacet"],
        #                                                     SpheNorm=False)
        #                                                     #,
        #                                                     #NFreqBands=self.VS.NFreqBands,
        #                                                     #DataCorrelationFormat=self.VS.StokesConverter.AvailableCorrelationProductsIds(),
        #                                                     #ExpectedOutputStokes=self.VS.StokesConverter.RequiredStokesProductsIds(),
        #                                                     #ListSemaphores=self.ListSemaphores)        

        SpheNorm = False
        FacetInfo = SM.DicoImager[iFacet]
        IDFacet=FacetInfo["IDFacet"]
        cf_dict=shared_dict.attach(SM.Path["cf_dict_path"])[IDFacet]
        #print iFacet,IDFacet
        GridMachine= ClassDDEGridMachine.ClassDDEGridMachine(SM.GD,
                                                             FacetInfo["DicoConfigGM"]["ChanFreq"],
                                                             FacetInfo["DicoConfigGM"]["NPix"],
                                                             lmShift=FacetInfo["lmShift"],
                                                             IDFacet=IDFacet,
                                                             SpheNorm=SpheNorm, 
                                                             NFreqBands=SM.NFreqBands,
                                                             DataCorrelationFormat=SM.AvailableCorrelationProductsIds,
                                                             ExpectedOutputStokes=SM.RequiredStokesProductsIds,
                                                             cf_dict=cf_dict)

        return GridMachine

    def InitGM(self,SM):
        self.DicoGM={}
        if SM.Type=="Hybrid":
            SM=SM.LSM[0]
        LFacets=SM.DicoImager.keys()
        for iFacet in LFacets:
            #print "Initialising Facet %i"%iFacet
            self.DicoGM[iFacet]=self.GiveGM(iFacet,SM)


    def predictKernelPolCluster(self,DicoData,SM,indKills=None,**kwargs):
        if self.SM is not None:
            SM=self.SM
        
        if indKills is not None and SM.Type!="Image":
            SM.SelectSubCat(indKills)
            
        if SM.Type=="Catalog":
            R=self.predictKernelPolClusterCatalog(DicoData,SM,**kwargs)
        elif SM.Type=="Image":
            R=self.predictKernelPolClusterImage(DicoData,SM,**kwargs)
        else:
            stop
            
        if indKills is not None and SM.Type!="Image":
            SM.RestoreCat()
            
        return R
            
    def predictKernelPolClusterCatalog(self,DicoData,SM,iDirection=None,ApplyJones=None,ApplyTimeJones=None,Noise=None,freq=None):
        if self.SM is not None:
            SM=self.SM
        self.DicoData=DicoData
        self.SourceCat=SM.SourceCat

        if freq is None:
            freq=DicoData["freqs_full"]
            
        times=DicoData["times"]
        nf=freq.size
        na=DicoData["infos"][0]
        
        nrows=DicoData["A0"].size
        DataOut=np.zeros((nrows,nf,4),self.CType)
        if nrows==0: return DataOut
        
        self.freqs=freq
        self.wave=299792458./self.freqs
        
        if iDirection!=None:
            ListDirection=[iDirection]
            #ListDirection=[SM.Dirs[iDirection]]
        else:
            ListDirection=SM.Dirs#range(SM.NDir)
        
        A0=DicoData["A0"]
        A1=DicoData["A1"]
        Map_VisToJones_Time=DicoData.get("Map_VisToJones_Time",None)
        Map_VisToJones_Freq=DicoData.get("Map_VisToJones_Freq",None)
        
        if ApplyJones is not None:
            #print "!!!!!",ApplyJones.shape
            #print "!!!!!",ApplyJones.shape
            #print "!!!!!",ApplyJones.shape
            na,NDir,_=ApplyJones.shape
            Jones=np.swapaxes(ApplyJones,0,1)
            Jones=Jones.reshape((NDir,na,4))
            JonesH=ModLinAlg.BatchH(Jones)

        TSmear=0.
        FSmear=0.
        DT=DicoData["infos"][1]
        UVW_dt=DicoData["uvw"]
        if self.DoSmearing:
            if "T" in self.DoSmearing:
                TSmear=1.
                UVW_dt=DicoData["UVW_dt"]
            if "F" in self.DoSmearing:
                FSmear=1.




        # self.SourceCat.m[:]=0
        # self.SourceCat.l[:]=0.1
        # self.SourceCat.I[:]=10
        # self.SourceCat.alpha[:]=0

        # DataOut=DataOut[1:2]
        # self.DicoData["uvw"]=self.DicoData["uvw"][1:2]
        # self.DicoData["A0"]=self.DicoData["A0"][1:2]
        # self.DicoData["A1"]=self.DicoData["A1"][1:2]
        # self.DicoData["IndexTimesThisChunk"]=self.DicoData["IndexTimesThisChunk"][1:2]
        # self.SourceCat=self.SourceCat[0:1]

        
        ColOutDir=np.zeros(DataOut.shape,np.complex64)

        for iCluster in ListDirection:
            ColOutDir.fill(0)

            indSources=np.where(self.SourceCat.Cluster==iCluster)[0]
            T=ClassTimeIt("  predict_Numexpr5")
            T.disable()
            ### new
            SourceCat=self.SourceCat[indSources].copy()
            #l=np.ones((1,),dtype=np.float64)#,float64(SourceCat.l).copy()
            l=np.require(SourceCat.l, dtype=np.float64, requirements=["A","C"])
            m=np.require(SourceCat.m, dtype=np.float64, requirements=["A","C"])
            
            #m=SourceCat.m#np.float64(SourceCat.m).copy()
            I=np.float32(SourceCat.I)
            Gmaj=np.float32(SourceCat.Gmaj)
            Gmin=np.float32(SourceCat.Gmin)
            GPA=np.float32(SourceCat.Gangle)
            alpha=np.float32(SourceCat.alpha)
            WaveL=np.float64(299792458./self.freqs)
            WaveL=np.require(WaveL, dtype=np.float64, requirements=["A","C"])

            flux=np.float32(SourceCat.I)
            alpha=SourceCat.alpha
            SType=SourceCat.Type
            
            #print("!!!!!!!!!!!!!!!!!!!",SType)
            
            dnu=np.float32(self.DicoData["dfreqs_full"])
            f0=(self.freqs/SourceCat.RefFreq[0])
            fluxFreq=np.float32(flux.reshape((flux.size,1))*(f0.reshape((1,f0.size)))**(alpha.reshape((alpha.size,1))))

            # print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            # d=15.*300/3600*np.pi/180
            # Gmin.fill(d)
            # Gmaj.fill(d)

            
            LSM=[l,m,fluxFreq,Gmin,Gmaj,GPA,SType]
            LFreqs=[WaveL,np.float32(self.freqs),dnu]

            LUVWSpeed=[UVW_dt,DT]

            LSmearMode=[FSmear,TSmear]
            #print("!!!!!!!!!!!",LSmearMode)
            #stop
            T.timeit("init")

            AllowEqualiseChan=IsChanEquidistant(DicoData["freqs_full"])

            if ApplyTimeJones is not None:
                # predict.predictJones(ColOutDir,(DicoData["uvw"]),LFreqs,LSM,LUVWSpeed,LSmearMode,ParamJonesList)

                # ColOutDir.fill(0)
                # predict.predictJones(ColOutDir,(DicoData["uvw"]),LFreqs,LSM,LUVWSpeed,LSmearMode,ParamJonesList,0)
                # d0=ColOutDir.copy()
                # ColOutDir.fill(0)

                # predict.predictJones(ColOutDir,(DicoData["uvw"]),LFreqs,LSM,LUVWSpeed,LSmearMode,ParamJonesList,AllowEqualiseChan)
                # ColOutDir0=ColOutDir.copy()
                # ColOutDir.fill(0)

                # predict.predictJones2(ColOutDir,(DicoData["uvw"]),LFreqs,LSM,LUVWSpeed,LSmearMode,ParamJonesList,AllowEqualiseChan)
                # print LSmearMode
                
                # print(LSmearMode)
                # print(SM.StrPhaseCenter_RADEC)
                predict.predictJones2_Gauss(ColOutDir,(DicoData["uvw"]),LFreqs,LSM,LUVWSpeed,LSmearMode,AllowEqualiseChan,
                                            self.LExp,
                                            self.LSinc,
                                            self.LBessel)
                
                T.timeit("predict")

                ParamJonesList=self.GiveParamJonesList(ApplyTimeJones,A0,A1,Map_VisToJones_Time,Map_VisToJones_Freq)
                J=ParamJonesList[-2]
                #print("JJJ",J[Map_VisToJones_Time,0,0,0,0,0][0])
                
                ra1,dec1=ApplyTimeJones["DicoFacetDirs"]["ra"],ApplyTimeJones["DicoFacetDirs"]["dec"]
                ra0,dec0=SM.ClusterCat.ra[iCluster],SM.ClusterCat.dec[iCluster]
                d=AngDist(ra0,ra1,dec0,dec1)
                iPreApplyDirs=np.argmin(d)
                #print("  Catalog:",iCluster,iPreApplyDirs,ra1.shape)
                JJ=ParamJonesList[-2]
                # print(JJ.shape,iPreApplyDirs)
                ParamJonesList=ParamJonesList+[iPreApplyDirs]
                
                predict.ApplyJones(ColOutDir,ParamJonesList)
                T.timeit("apply")


                # print ColOutDir

                #d1=ColOutDir.copy()
                #ind=np.where(d0!=0)
                #print np.max((d0-d1)[ind]/(d0[ind]))
                #stop
            else:
                #predict.predict(ColOutDir,(DicoData["uvw"]),LFreqs,LSM,LUVWSpeed,LSmearMode)
                #AllowEqualiseChan=0
                #predict.predict(ColOutDir,(DicoData["uvw"]),LFreqs,LSM,LUVWSpeed,LSmearMode,AllowEqualiseChan)
                #d0=ColOutDir.copy()
                #ColOutDir.fill(0)


                predict.predictJones2_Gauss(ColOutDir,(DicoData["uvw"]),LFreqs,LSM,LUVWSpeed,LSmearMode,AllowEqualiseChan,
                                            self.LExp,
                                            self.LSinc,
                                            self.LBessel)
                #print ColOutDir
                #predict.predict(ColOutDir,(DicoData["uvw"]),LFreqs,LSM,LUVWSpeed,LSmearMode,AllowEqualiseChan)
                T.timeit("predict")
                # d0=ColOutDir.copy()
                # ColOutDir.fill(0)

                # predict_np19.predict(ColOutDir,(DicoData["uvw"]),LFreqs,LSM,LUVWSpeed,LSmearMode,AllowEqualiseChan)
                # print ColOutDir,d0
                # d1=ColOutDir.copy()
                # ind=np.where(d0!=0)
                # print np.max((d0-d1)[ind]/(d0[ind]))
                # stop


            del(l,m,I,SourceCat,alpha,WaveL,flux,dnu,f0,fluxFreq,LSM,LFreqs)
            T.timeit("del")


                # d1=ColOutDir
                # ind=np.where(d0!=0)
                # print np.max((d0-d1)[ind]/(d0[ind]))
                # stop
                



            if Noise!=None:
                ColOutDir+=Noise*(np.random.randn(*ColOutDir.shape)+1j*np.random.randn(*ColOutDir.shape))
                stop
            DataOut+=ColOutDir
            T.timeit("add")
        #del(LFreqs,LSM,LUVWSpeed,LSmearMode)
        #del(ColOutDir)

        return DataOut


    ######################################################
    ######################################################
    ######################################################


    def predictKernelPolClusterImage(self,DicoData,SM,iDirection=None,ApplyJones=None,ApplyTimeJones=None,Noise=None,ForceNoDecorr=False):
        if self.SM is not None:
            SM=self.SM

        T=ClassTimeIt("predictKernelPolClusterImage")
        T.disable()
        self.DicoData=DicoData

        
        
        freq=DicoData["freqs_full"]
        times=DicoData["times"]
        nf=freq.size
        na=DicoData["infos"][0]
        
        nrows=DicoData["A0"].size
        T.timeit("0")
        DataOut=np.zeros((nrows,nf,4),self.CType)
        if nrows==0: return DataOut
        
        self.freqs=freq
        self.wave=299792458./self.freqs
        
        if iDirection!=None:
            ListDirection=[iDirection]
        else:
            ListDirection=SM.Dirs#range(SM.NDir)
        
        A0=DicoData["A0"]
        A1=DicoData["A1"]
        Map_VisToJones_Time=DicoData.get("Map_VisToJones_Time",None)
        Map_VisToJones_Freq=DicoData.get("Map_VisToJones_Freq",None)
        # if ApplyJones!=None:
        #     na,NDir,_=ApplyJones.shape
        #     Jones=np.swapaxes(ApplyJones,0,1)
        #     Jones=Jones.reshape((NDir,na,4))
        #     JonesH=ModLinAlg.BatchH(Jones)

        TSmear=0.
        FSmear=0.


        if self.DoSmearing!=0:
            if "T" in self.DoSmearing:
                TSmear=1.
            if "F" in self.DoSmearing:
                FSmear=1.

        # self.SourceCat.m[:]=0
        # self.SourceCat.l[:]=0.1
        # self.SourceCat.I[:]=10
        # self.SourceCat.alpha[:]=0

        # DataOut=DataOut[1:2]
        # self.DicoData["uvw"]=self.DicoData["uvw"][1:2]
        # self.DicoData["A0"]=self.DicoData["A0"][1:2]
        # self.DicoData["A1"]=self.DicoData["A1"][1:2]
        # self.DicoData["IndexTimesThisChunk"]=self.DicoData["IndexTimesThisChunk"][1:2]
        # self.SourceCat=self.SourceCat[0:1]

        DT=DicoData["infos"][1]
        #UVW_dt=DicoData["UVW_dt"]
        Dnu=DicoData["dfreqs_full"][0]

        ColOutDir=np.zeros(DataOut.shape,np.complex64)
        F_zeroPredict=np.zeros(DataOut.shape,bool)
        DATA=DicoData

        T.timeit("1")
        ListFacets=SM.givePreApplyDirs(iDirection)

        for iFacet in ListFacets:
            if SM.DicoImager[iFacet]["SumFlux"]==0: continue
            GridMachine=self.DicoGM[iFacet]#self.GiveGM(iFacet,SM)
            T.timeit("2: GM")
            uvwThis=DATA["uvw"]
            
            flagsThis=DATA.get("flags_image",np.zeros(DataOut.shape,np.bool8))
            times=DATA["times"]
            A0=DATA["A0"]
            A1=DATA["A1"]
            A0A1=A0,A1
            freqs=DATA["freqs_full"]


            DicoJonesMatrices=None
            #DicoJonesMatrices=ApplyTimeJones



            # ModelSharedMemName="%sModelImage.Facet_%3.3i"%(self.IdSharedMem,iFacet)
            # print "Facet %i: take model image %s"%(iFacet,ModelSharedMemName)
            # ModelIm = NpShared.GiveArray(ModelSharedMemName)

            #ModelIm = NpShared.UnPackListArray("%sGrids"%self.IdSharedMem)[iFacet]
            ModelIm = SM._model_dict[iFacet]["FacetGrid"]
            
            ChanMapping=np.int32(SM.ChanMappingDegrid)
            # print ChanMapping
            
            GridMachine.LSmear=[]
            DecorrMode = SM.GD["RIME"]["DecorrMode"]
            CondSmear=(not ForceNoDecorr) and (('F' in DecorrMode) | ("T" in DecorrMode))
            if CondSmear:
                #print "DOSMEAR",ForceNoDecorr, (('F' in DecorrMode) | ("T" in DecorrMode))
                uvw_dt = DicoData["UVW_dt"]#DATA["uvw_dt"]
                lm_min=None
                if SM.GD["RIME"]["DecorrLocation"]=="Edge":
                    lm_min=SM.DicoImager[iFacet]["lm_min"]
                    
                #lm_PhaseCenter=DATA["lm_PhaseCenterflags_image"]
                
                GridMachine.setDecorr(uvw_dt, DT, Dnu, SmearMode=SM.GD["RIME"]["DecorrMode"],
                                      lm_min=lm_min,
                                      lm_PhaseCenter=SM.lm_PhaseCenter)
                

            T.timeit("2: Stuff")

            # print """
            # print times.shape,times.dtype
            # print uvwThis.shape,uvwThis.dtype
            # print ColOutDir.shape,ColOutDir.dtype
            # print flagsThis.shape,flagsThis.dtype
            # print ModelIm.shape,ModelIm.dtype"""
            
            # print times.shape,times.dtype
            # print uvwThis.shape,uvwThis.dtype
            # print ColOutDir.shape,ColOutDir.dtype
            # print flagsThis.shape,flagsThis.dtype
            # print ModelIm.shape,ModelIm.dtype
            # stop
            vis=GridMachine.get(times,uvwThis,ColOutDir,flagsThis,A0A1,ModelIm,DicoJonesMatrices=DicoJonesMatrices,freqs=freqs,
                                ImToGrid=False,ChanMapping=ChanMapping)

            T.timeit("2: Predict")
            # get() is substracting
            if ApplyTimeJones is not None and self._BeamAtFacet:

                #print("  Image: Facet",iDirection,iFacet)
                ParamJonesList=self.GiveParamJonesList(ApplyTimeJones,A0,A1,Map_VisToJones_Time,Map_VisToJones_Freq)
                ParamJonesList=ParamJonesList+[iFacet]

                # print "facet"
                # import killMS.Other.rad2hmsdms
                # RA,DEC=SM.DicoImager[iFacet]["RaDec"]
                # sra=killMS.Other.rad2hmsdms.rad2hmsdms(RA,Type="ra")
                # sdec=killMS.Other.rad2hmsdms.rad2hmsdms(DEC,Type="dec")
                # print iFacet,sra,sdec
                
                predict.ApplyJones(ColOutDir,ParamJonesList)

            F_zeroPredict[ColOutDir==0]=1
            DataOut-=ColOutDir
            ColOutDir.fill(0)
            
        if ApplyTimeJones is not None and not self._BeamAtFacet:
            #print "tessel"
            #print "apply in direction %i"%iDirection
            #print("  Image: ApplyTessel",iDirection)
            ParamJonesList=self.GiveParamJonesList(ApplyTimeJones,A0,A1,Map_VisToJones_Time,Map_VisToJones_Freq)
            ParamJonesList=ParamJonesList+[iDirection]
            predict.ApplyJones(DataOut,ParamJonesList)
            T.timeit("apply")
            
        if np.count_nonzero(np.isnan(DataOut)):
            stop

        T.timeit("2: End")

        DataOut[F_zeroPredict]=0
        return DataOut


#####################################################


class ClassPredictParallel():
    def __init__(self,Precision="S",NCPU=6,IdMemShared="",DoSmearing=False,BeamAtFacet=False):
        self.NCPU=NCPU
        ne.set_num_threads(self.NCPU)
        if Precision=="D":
            self.CType=np.complex128
            self.FType=np.float64
        if Precision=="S":
            self.CType=np.complex64
            self.FType=np.float32
        self.IdMemShared=IdMemShared
        self.DoSmearing=DoSmearing
        self._BeamAtFacet = BeamAtFacet
        self.SM=None
        self.PM=ClassPredict(NCPU=1,DoSmearing=DoSmearing,BeamAtFacet=BeamAtFacet)



    def setSM(self,SM):
        self.SM=SM
        self.PM.setSM(SM)

    def InitGM(self,SM):
        if SM.Type=="Image":
            self.PM.InitGM(SM)
        elif SM.Type=="Hybrid":
            self.PM.InitGM(SM.LSM[0])
            
        
    def GiveCovariance(self,DicoDataIn,ApplyTimeJones,SM,Mode="DDECovariance"):
        
        if not isinstance(DicoDataIn,shared_dict.SharedDict):
            DicoData=dict_to_shm("DicoMemChunk",DicoDataIn)
        else:
            DicoData=DicoDataIn
            
            
        NameApplyTimeJones=None
        if ApplyTimeJones is not None:
            if not isinstance(ApplyTimeJones,shared_dict.SharedDict):
                ApplyTimeJones=shared_dict.dict_to_shm("ApplyTimeJones",ApplyTimeJones)
            NameApplyTimeJones=ApplyTimeJones.path

        nrow,nch,_=DicoData["data"].shape
        
        RowList=np.int64(np.linspace(0,nrow,self.NCPU+1))
        row0=RowList[0:-1]
        row1=RowList[1::]
        
        NJobs=row0.size
        for iJob in range(NJobs):
            NameJob="workerPredict:%i:%i"%(row0[iJob],row1[iJob])
            APP.APP.runJob(NameJob,
                           self._workerPredict,
                           args=(None,
                                 Mode,
                                 0,
                                 self.DoSmearing,
                                 self._BeamAtFacet,
                                 row0[iJob],
                                 row1[iJob],
                                 None,
                                 DicoData.path,
                                 NameApplyTimeJones))#,serial=True)
            
        Message=APP.APP.awaitJobResults("workerPredict*",progress="Compute Weights")

        DicoDataIn["W"]=DicoData["W"]


    def ApplyCal(self,DicoDataIn,ApplyTimeJones,iCluster):

        if not isinstance(DicoDataIn,shared_dict.SharedDict):
            DicoData=dict_to_shm("DicoMemChunk",DicoDataIn)
        else:
            DicoData=DicoDataIn

        NameApplyTimeJones=None
        if ApplyTimeJones is not None:
            if not isinstance(ApplyTimeJones,shared_dict.SharedDict):
                ApplyTimeJones=shared_dict.dict_to_shm("ApplyTimeJones",ApplyTimeJones)
            NameApplyTimeJones=ApplyTimeJones.path
            
        nrow,nch,_=DicoData["data"].shape
        
        RowList=np.int64(np.linspace(0,nrow,self.NCPU+1))
        row0=RowList[0:-1]
        row1=RowList[1::]
        
        NJobs=row0.size
        Mode="ApplyCal"
        for iJob in range(NJobs):
            NameJob="workerPredict:%i:%i"%(row0[iJob],row1[iJob])
            APP.APP.runJob(NameJob,
                           self._workerPredict,
                           args=(None,
                                 Mode,
                                 iCluster,
                                 self.DoSmearing,
                                 self._BeamAtFacet,
                                 row0[iJob],
                                 row1[iJob],
                                 None,
                                 DicoData.path,
                                 NameApplyTimeJones))#,serial=True)

        Message=APP.APP.awaitJobResults("workerPredict*",progress="Apply solutions")

        DicoDataIn["data"]=DicoData["data"]


    def predictKernelPolCluster(self,DicoData,SM,iDirection=None,ApplyJones=None,ApplyTimeJones=None,Noise=None,
                                indKills=None,Progress=None):

        
        T=ClassTimeIt("predict")
        T.disable()
        
        if not isinstance(DicoData,shared_dict.SharedDict):
            DicoData=dict_to_shm("DicoMemChunk",DicoData)
            
        T.timeit("attach shared dict 0")
        NameApplyTimeJones=None
        if ApplyTimeJones is not None:
            if not isinstance(ApplyTimeJones,shared_dict.SharedDict):
                ApplyTimeJones=shared_dict.dict_to_shm("ApplyTimeJones",ApplyTimeJones)
            NameApplyTimeJones=ApplyTimeJones.path
            
        T.timeit("attach shared dict 1")
        
        nrow,nch,_=DicoData["data"].shape
        
        RowList=np.int64(np.linspace(0,nrow,(self.NCPU+1)*5))
        row0=RowList[0:-1]
        row1=RowList[1::]
        
        T.timeit("linspace")

        nr,nch,npol=DicoData["data"].shape
        #shared_dict.delDict("DicoPredictData")
        #DicoPredictData=shared_dict.create("DicoPredictData")

        if iDirection is None:
            LDirections=[None]
        else:
            if not isinstance(iDirection,list):
                LDirections=[iDirection]
            else:
                LDirections=iDirection
                if "iDirMap" not in DicoData.keys():
                    DicoData.addSubdict("iDirMap")
                for iDir,iDirIndex in enumerate(LDirections):
                    DicoData["iDirMap"][iDirIndex]=iDir

        if "PredictArray" not in DicoData.keys():
            #print("kjhdfkljjkljklh")
            #print("kjhdfkljjkljklh")
            #print("kjhdfkljjkljklh")
            DicoData["PredictArray"]=np.zeros((len(LDirections),nr,nch,npol),DicoData["data"].dtype)
            
        # DicoData["PredictArray"].fill(0)
        
        T.timeit("DicoPredictData0")
            
        #PredictArray=DicoPredictData["PredictArray"]
        #PredictArray.fill(0)
        
        # T.timeit("DicoPredictData1")
        
        NJobs=row0.size
        Mode="Predict"
        T.timeit("Init1")
        for iJob in range(NJobs):
            for iDirection in LDirections:
                iDirTag=0
                if iDirection is not None: iDirTag=iDirection

                NameJob="workerPredict:%i:%i%i"%(row0[iJob],row1[iJob],iDirTag)
                APP.APP.runJob(NameJob,
                               self._workerPredict,
                               args=(None,
                                     Mode,
                                     iDirection,
                                     self.DoSmearing,
                                     self._BeamAtFacet,
                                     row0[iJob],
                                     row1[iJob],
                                     indKills,
                                     DicoData.path,
                                     NameApplyTimeJones))#,serial=True)

        T.timeit("runJobs")
        if not Progress:
            Message=APP.APP.awaitJobResults("workerPredict*")
        else:
            Message=APP.APP.awaitJobResults("workerPredict*",progress="Predict")
        T.timeit("await ok")

        if len(LDirections)==1:
            PredictArray=DicoData["PredictArray"][0]
        else:
            PredictArray=DicoData["PredictArray"]
        T.timeit("await ok")
            
        return PredictArray



    def _workerPredict(self,SM=None,
                       Mode="Predict",
                       iDirection=None,
                       DoSmearing=False,
                       BeamAtFacet=False,
                       Row0=None,
                       Row1=None,
                       indKills=None,
                       NameDicoData=None,
                       NameApplyTimeJones=None):
        
        D=shared_dict.attach(NameDicoData)
        #print(NameDicoData)
        
        DicoData={}
        DicoData["data"]=D["data"][Row0:Row1]
        DicoData["flags"]=D["flags"][Row0:Row1]
        DicoData["A0"]=D["A0"][Row0:Row1]
        DicoData["A1"]=D["A1"][Row0:Row1]
        DicoData["times"]=D["times"][Row0:Row1]
        DicoData["uvw"]=D["uvw"][Row0:Row1]
        DicoData["freqs"]=D["freqs"]
        DicoData["dfreqs"]=D["dfreqs"]
        DicoData["freqs_full"]=D["freqs_full"]
        DicoData["dfreqs_full"]=D["dfreqs_full"]
        # DicoData["UVW_dt"]=D["UVW_dt"]
        DicoData["infos"]=D["infos"]

        if "Map_VisToJones_Time" in D.keys():
            DicoData["Map_VisToJones_Time"]=D["Map_VisToJones_Time"][Row0:Row1]
            DicoData["Map_VisToJones_Freq"]=D["Map_VisToJones_Freq"]
        
        #print(DicoData["Map_VisToJones_Time"])
        #DicoData["IndRows_All_UVW_dt"]=D["IndRows_All_UVW_dt"]
        #DicoData["All_UVW_dt"]=D["All_UVW_dt"]
        if self.DoSmearing and "T" in self.DoSmearing:
            DicoData["UVW_dt"]=D["UVW_dt"][Row0:Row1]

        # DicoData["IndexTimesThisChunk"]=D["IndexTimesThisChunk"][Row0:Row1]
        # it0=np.min(DicoData["IndexTimesThisChunk"])
        # it1=np.max(DicoData["IndexTimesThisChunk"])+1
        # DicoData["UVW_RefAnt"]=D["UVW_RefAnt"][it0:it1,:,:]

        if "W" in D.keys():
            DicoData["W"]=D["W"][Row0:Row1]

        if "resid" in D.keys():
            DicoData["resid"]=D["resid"][Row0:Row1]

        # ApplyTimeJones_shared=shared_dict.attach("ApplyTimeJones")
        # ApplyTimeJones={}
        # for k in ApplyTimeJones_shared.keys():
        #     if k=="Map_VisToJones_Time":
        #         pass#ApplyTimeJones[k]=ApplyTimeJones_shared[k].copy()
        #         #stop
        #     else:
        #         ApplyTimeJones[k]=ApplyTimeJones_shared[k]
        
        ApplyTimeJones=None
        if NameApplyTimeJones is not None:
            ApplyTimeJones=shared_dict.attach(NameApplyTimeJones)
        #print(NameApplyTimeJones)
        #JonesMatrices=ApplyTimeJones["Beam"]
        #print ApplyTimeJones["Beam"].flags
        
            
        PM=self.PM#ClassPredict(NCPU=1,DoSmearing=self.DoSmearing,BeamAtFacet=self._BeamAtFacet)
        
        #print DicoData.keys()


        if Mode=="Predict":
            # if Row0==17980: stop
            PredictData=PM.predictKernelPolCluster(DicoData,
                                                   SM,
                                                   iDirection=iDirection,
                                                   ApplyTimeJones=ApplyTimeJones,
                                                   indKills=indKills)
            
            
            if iDirection is None:
                iDir=0
            else:
                iDir=D["iDirMap"][iDirection]
            PredictArray=D["PredictArray"][iDir]
            PredictArray[Row0:Row1]=PredictData[:]
        elif Mode=="ApplyCal":
            PM.ApplyCal(DicoData,ApplyTimeJones,iDirection)
        elif Mode=="DDECovariance":
            PM.GiveCovariance(DicoData,ApplyTimeJones,SM)
        elif Mode=="ResidAntCovariance":
            PM.GiveResidAntCovariance(DicoData,ApplyTimeJones,SM)

        return True
