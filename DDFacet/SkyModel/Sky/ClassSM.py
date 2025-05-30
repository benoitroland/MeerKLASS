from __future__ import division, absolute_import, print_function


import numpy as np
from . import ModTigger
from . import ModSMFromNp
from . import ModSMFromFITS
from DDFacet.Other import logger
log=logger.getLogger("ClassSM")

from SkyModel.Other import rad2hmsdms
from SkyModel.Other import ModColor
from SkyModel.Array import RecArrayOps
from SkyModel.Sky.ClassClusterClean import ClassClusterClean
from SkyModel.Sky import ClassAppendSource
from SkyModel.Sky.ClassClusterTessel import ClassClusterTessel
from SkyModel.Sky.ClassClusterRadial import ClassClusterRadial
from SkyModel.Sky.ClassClusterSquareRadial import ClassClusterSquareRadial
from SkyModel.Sky.ClassClusterKMean import ClassClusterKMean
#import ModClusterRadial
from pyrap.images import image
import scipy.linalg
from SkyModel.Sky.ModBBS2np import ReadBBSModel
#from SkyModel.Sky.ModSMFromFITS import ReadSMFromFITS
import DDFacet.Other.MyPickle
from DDFacet.ToolsDir.ModToolBox import EstimateNpix
from DDFacet.ToolsDir.rad2hmsdms import rad2hmsdms

from . import ModRegFile
import time
from DDFacet.ToolsDir import ModCoord

dtypeSourceList=[('Name','|S200'),('Patch','|S200'),('ra',np.float32),('dec',np.float32),('Sref',np.float32),('I',np.float32),('Q',np.float32),\
                 ('U',np.float32),('V',np.float32),('RefFreq',np.float32),('alpha',np.float32),('ESref',np.float32),\
                 ('Ealpha',np.float32),('kill',np.int32),('Cluster',np.int32),('Type',np.int32),('Gmin',np.float32),\
                 ('Gmaj',np.float32),('Gangle',np.float32),("Select",np.int32),('l',np.float32),('m',np.float32),("Exclude",bool)]

dtypeClusterCat=[('Name','|S200'),('ra',np.float32),('dec',np.float32),('SumI',np.float32),("Cluster",int),('l',np.float32),('m',np.float32)]

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


class ClassSMConcat():
    def __init__(self,LSM):
        self.Type="Hybrid"
        self.LSM=LSM


    def ComputeMapping(self):
        self.NDir=np.sum([SM.NDir for SM in self.LSM])
        self.NDirMain=self.NDir
        self.ClusterCat=np.concatenate([SM.ClusterCat for SM in self.LSM])
        self.ClusterCat=self.ClusterCat.view(np.recarray)
        self.iDir_to_SMiDir={}
        iDir=0
        for SM in self.LSM:
            for iDirThisSM in range(SM.NDir):
                self.iDir_to_SMiDir[iDir]=(SM,iDirThisSM)
                iDir+=1
        self.MapClusterCatOrigToCut=np.concatenate([SM.MapClusterCatOrigToCut for SM in self.LSM])
        self.NDirsOrig=np.sum([SM.NDirsOrig for SM in self.LSM])
        self.ClusterCatOrig=np.concatenate([SM.ClusterCatOrig for SM in self.LSM])
        self.ClusterCatOrig=self.ClusterCatOrig.view(np.recarray)
        log.print("Parameter space has %i directions (of %i before flux cut)"%(self.NDir,self.NDirsOrig))

    # def setDicoJonesDirToPreApplyDirs(self,radec):
    #     for SM in self.LSM:
    #         SM.setDicoJonesDirToPreApplyDirs(radec)


        
    def Calc_LM(self,rac,decc):
        for SM in self.LSM:
            if SM.Type=="Catalog":
                SM.Calc_LM(rac,decc)

    def give_SM_iDir(self,iDir):
        return self.iDir_to_SMiDir[iDir]



class ClassSM():
    def __init__(self,infile,infile_cluster="",killdirs=[],invert=False,DoPrintCat=False,\
                 ReName=False,
                 DoREG=False,SaveNp=False,NCluster=0,DoPlot=True,Tigger=False,\
                 FromExt=None,ClusterMethod=1,SelSource=False,DoPrint=True,ListBody=None,
                 radecCenter=None):
        self.ClusterCatExt=None
        self.ClusterMethod=ClusterMethod
        self.infile_cluster=infile_cluster
        self.TargetList=infile
        self.Type="Catalog"
        self.DoPrint=DoPrint
        self.D_FITS=None
        self.DDF_GD=None
        self.InputCatIsEmpty=False
        self.ClusterCat=None
        self.ModelIsClustered=False
        if (type(infile).__name__=="instance") or (type(infile).__name__=="ClassImageSM") or (type(infile).__name__=="ClassSM"):
            Cat=infile.SourceCat.copy()
            Cat=Cat.view(np.recarray)
            self.DoPrint=0
        elif infile.endswith(".npy"):
            Cat=np.load(infile)
            Cat=Cat.view(np.recarray)
        elif infile=="Empty":
            Cat=np.zeros((0,),dtype=dtypeSourceList)
            Cat=Cat.view(np.recarray)
            self.InputCatIsEmpty=True
        elif infile.endswith(".pickle"):
            D=DDFacet.Other.MyPickle.Load(infile)
            Cat0=D["SourceCat"].view(np.recarray)
            self.NSources=Cat0.shape[0]
            Cat=np.zeros((self.NSources,),dtype=dtypeSourceList)
            for field in list(Cat0.dtype.fields.keys()):
                Cat[field][:]=Cat0[field][:]
            Cat=Cat.view(np.recarray)
            self.ClusterCat=D["ClusterCat"]
            self.D_FITS=D["D_FITS"]
            self.DDF_GD=D["DDF_GD"]
        elif Tigger:
            Cat=ModTigger.ReadTiggerModel(infile)
        elif FromExt is not None:
            Cat=ModSMFromNp.ReadFromNp(FromExt)
        elif ".fits"==infile[-5:]:
            Cat,self.D_FITS=ModSMFromFITS.ReadSMFromFITS(infile)
        else:
            Cat,self.ModelIsClustered=ReadBBSModel(infile,infile_cluster=infile_cluster)
            
            Cat.Gmaj[Cat.Type==2]*=2.*np.sqrt(2.*np.log(2))
            Cat.Gmin[Cat.Type==2]*=2.*np.sqrt(2.*np.log(2))

        self.DicoJonesDirToPreApplyDirs=None
        self.SourceCat=Cat
        self.killdirs=killdirs
        self.invert=invert
        self.infile=infile
        self.REGFile=None
        if not self.InputCatIsEmpty:
            if self.ClusterCat is not None:
                self.Dirs=sorted(self.ClusterCat.tolist())
                self.NDir=len(self.Dirs)
            else:
                self.Dirs=sorted(list(set(self.SourceCat.Cluster.tolist())))
                self.NDir=np.max(self.SourceCat.Cluster)+1
            self.NSources=Cat.shape[0]
        else:
            self.Dirs=[]
            self.NDir=0
            self.NSources=0
            
        self.StrPhaseCenter_RADEC=self.rarad=self.decrad=None

        try:
            if radecCenter is not None:
                log.print(("Computing lm of the center at externaly specified rac/dec..."))
                self.rarad,self.decrad=radecCenter
            elif self.D_FITS is not None:
                print(("Reading lm of the center from FITS"))
                #self.rarad=SourceCat.ra[indIMax]#
                self.rarad=self.D_FITS["rac"]
                self.decrad=self.D_FITS["decc"]
            else:
                SourceCat=self.SourceCat
                indIMax=np.argmax(SourceCat.I)
                log.print(("Computing lm of the center at flux-weighted mean rac/dec"))
                #self.rarad=SourceCat.ra[indIMax]#
                self.rarad=np.sum(SourceCat.I*SourceCat.ra)/np.sum(SourceCat.I)
                self.decrad=np.sum(SourceCat.I*SourceCat.dec)/np.sum(SourceCat.I)
                
            self.CoordMachine = ModCoord.ClassCoordConv(self.rarad, self.decrad)
        except:
            pass

        if self.rarad is not None:
            self.StrPhaseCenter_RADEC  = rad2hmsdms(self.rarad,Type="ra").replace(" ",":") , rad2hmsdms(self.decrad,Type="dec").replace(" ",".")
            
        self.NDirMain=self.NDir
        
        if ListBody is not None:
            log.print("Append sources to the sky model:")
            CAS=ClassAppendSource.ClassAppendSource(self,ListBody)
            CAS.appendAll()
            self.ClusterCat=None

        if self.ClusterCat is None:
            self.Dirs=sorted(list(set(self.SourceCat.Cluster.tolist())))
            self.NDir=0
            if self.SourceCat.size>0:
                self.NDir=np.max(self.SourceCat.Cluster)+1
            self.NSources=Cat.shape[0]
            self.BuildClusterCat()
            self.Dirs=sorted(list(set(self.SourceCat.Cluster.tolist())))
            if self.SourceCat.size>0:
                self.NDir=np.max(self.SourceCat.Cluster)+1
            self.NSources=Cat.shape[0]


        self.SetSelection()

        if self.DoPrint: self.PrintBasics()


        #self.print_sm2()

    # def setDicoJonesDirToPreApplyDirs(self,radecPreapply):
    #     ra1,dec1=radecPreapply
    #     ra1,dec1=ra1.reshape((1,-1)),dec1.reshape((1,-1))
    #     ra0,dec0=self.ClusterCat.ra.reshape((-1,1)),self.ClusterCat.dec.reshape((-1,1))
    #     d=AngDist(ra0,ra1,dec0,dec1)
    #     iDirPreApply=np.argmin(d,axis=1)
    #     self.DicoJonesDirToPreApplyDirs=iDirPreApply
        
    # def givePreApplyDirs(self,iClusterDir):
        
    #     return self.DicoJonesDirToPreApplyDirs[iClusterDir]
        
    def SetSelection(self):
        Cat=self.SourceCat
        killdirs=self.killdirs
        invert=self.invert
        self.SourceCatKeepForSelector=self.SourceCat.copy()

        # if DoPrintCat:
        #     self.print_sm(Cat)

        self.SourceCat.kill=1
        if killdirs!=[]:
            self.SourceCat.kill=0
            for StrPiece in killdirs:
                if StrPiece=="":
                    continue
                elif "x" in StrPiece:
                    if np.all(self.SourceCat.l==0):
                        self.Calc_LM(self.rarad,self.decrad)
                    snpix,scell =StrPiece.split("x")
                    Npix, _ = EstimateNpix(float(snpix), Padding=1)
                    Cell=float(scell)/3600*np.pi/180
                    Cx0 = (self.SourceCat.l>(-Npix//2*Cell))
                    Cx1 = (self.SourceCat.l<(Npix//2*Cell))
                    Cy0 = (self.SourceCat.m>(-Npix//2*Cell))
                    Cy1 = (self.SourceCat.m<(Npix//2*Cell))
                    ind=np.where(Cx0&Cx1&Cy0&Cy1)[0]
                    self.SourceCat.kill[ind]=1
                else:
                    for i in range(len(self.SourceCat)):
                        Name=self.SourceCat.Name[i]
                        if "byte" in type(Name).__name__: Name=Name.decode("utf-8")
                        if StrPiece in Name: self.SourceCat.kill[i]=1
                    
                    
        if invert:
            ind0=np.where(self.SourceCat.kill==0)[0]
            ind1=np.where(self.SourceCat.kill==1)[0]
            self.SourceCat.kill[ind0]=1
            self.SourceCat.kill[ind1]=0

            
        # print(self.SourceCat.Name)
        # print(self.SourceCat.kill)
            
        #log.print(ModColor.Str())

    def Save(self):
        infile=self.infile
        print(ModColor.Str(" SkyModel PROPERTIES: "))
        print("   - SkyModel File Name: %s"%ModColor.Str(infile,col="green"))
        if self.REGFile is not None: print("   - ds9 region file: %s"%ModColor.Str(self.REGFile,col="green"))
        npext=""
        if not(".npy" in infile): npext=".npy"
        self.NpFile="%s%s"%(infile,npext)
        np.save(infile,self.SourceCat)

        FileClusterCat="%s.ClusterCat.npy"%(self.infile)
        print("   - ClusterCat File Name: %s"%ModColor.Str(FileClusterCat,col="green"))
        np.save(FileClusterCat,self.ClusterCat)

        self.PrintBasics()

    def SavePickle(self):
        infile=self.infile
        print(ModColor.Str(" SkyModel PROPERTIES: "))
        npext=""
        if not(".pickle" in infile): npext=".pickle"
        self.NpFile="%s%s"%(infile,npext)
        D={"SourceCat":self.SourceCat,
           "D_FITS":self.D_FITS,
           "DDF_GD":self.DDF_GD,
           "ClusterCat":self.ClusterCat}
        print("   - SkyModel File Name: %s"%ModColor.Str(self.NpFile,col="green"))
        DDFacet.Other.MyPickle.Save(D,self.NpFile)
        if self.REGFile is not None: print("   - ds9 region file: %s"%ModColor.Str(self.REGFile,col="green"))


        self.PrintBasics()

        
    def PrintBasics(self):
        infile=self.infile
        # if "instance" in str(type(infile)): return
        npext=""
        if not(".npy" in infile): npext=".npy"
        print("   - Numpy catalog file: %s"%ModColor.Str("%s%s"%(infile,npext),col="green"))

        #print "Oufile: %s"%self.infile_cluster
        #if infile_cluster!="":
        #print "   - Cluster File Name: %s"%self.infile_cluster
        print("   - Number of Sources  = ",self.SourceCat.shape[0])
        print("   - Number of Directions  = ",self.NDir)
        Np=np.count_nonzero(self.SourceCat.Type==0)
        Ng=np.count_nonzero(self.SourceCat.Type==1)
        Nb=np.count_nonzero(self.SourceCat.Type==2)
        print("   - Number of [ POINT | GAUSSIANS | BOX ] : [ %i | %i | %i ]"%(Np,Ng,Nb))
        NSub=np.count_nonzero(self.SourceCat.kill)
        print("   - Selected for subtraction: %i/%i"%(NSub,self.SourceCat.kill.size))
        print()

    def Cluster(self,NCluster=1,DoPlot=True,PreCluster="",FromClusterCat=""):

        if PreCluster!="":
            R=ModRegFile.RegToNp(PreCluster)
            R.Read()
            #R.Cluster()
            PreClusterCat=R.CatSel
            ExcludeCat=R.CatExclude
        else:
            PreClusterCat=None
            ExcludeCat=None

        if ExcludeCat is not None:
            for j in range(ExcludeCat.shape[0]):
                d=np.sqrt((self.SourceCat.ra-ExcludeCat.ra[j])**2+(self.SourceCat.dec-ExcludeCat.dec[j])**2)
                self.SourceCat.Exclude[d<ExcludeCat.Radius[j]]=True

        if FromClusterCat and (NCluster==0):
            self.ClusterCatExt=self.ClusterCat=np.load(FromClusterCat).view(np.recarray)
            
            NCluster=self.ClusterCat.shape[0]
            RA=self.SourceCat.ra
            DEC=self.SourceCat.dec
            RAd=self.ClusterCat.ra
            DECd=self.ClusterCat.dec
            
            d=AngDist(RA.reshape((-1,1)),RAd.reshape((1,-1)),
                      DEC.reshape((-1,1)),DECd.reshape((1,-1)))
            C=np.argmin(d,axis=1)
            self.SourceCat.Cluster[:]=C[:]
            self.Rename()
            self.Dirs=sorted(list(set(self.SourceCat.Cluster.tolist())))
            self.NDir=self.ClusterCat.shape[0]
            for iDir in self.Dirs:
                ind=np.where(self.SourceCat.Cluster==iDir)[0]
                cat=self.SourceCat[ind]
                self.ClusterCat.SumI[iDir]=np.sum(cat.I)
            NZeroFlux=np.count_nonzero(self.ClusterCat.SumI==0)
            if NZeroFlux!=0:
                log.print("There are %i zero flux directions (on %i directions)"%(NZeroFlux,self.NDir))
                ind=np.where(self.ClusterCat.SumI==0)[0]
                for iDirNull in ind:
                    c=self.SourceCat[0:1].copy()
                    c=c.view(np.recarray)
                    c.ra[0]=self.ClusterCat.ra[iDirNull]
                    c.dec[0]=self.ClusterCat.ra[iDirNull]
                    c.Sref[0]=0
                    c.I[0]=0
                    c.Gmaj[0]=0
                    c.Gmin[0]=0
                    c.Type[0]=0
                    c.Cluster[0]=iDirNull
                    self.SourceCat=np.concatenate([self.SourceCat,c],axis=0)
                    self.SourceCat=self.SourceCat.view(np.recarray)

            return
                
        if NCluster==0:
            self.SourceCat.Cluster=np.arange(self.SourceCat.shape[0])
        elif NCluster==1:
            self.SourceCat.Cluster=0
        else:
            self.cluster(NCluster,DoPlot,PreClusterCat=PreClusterCat,FromClusterCat=FromClusterCat)#,ExcludeCat=ExcludeCat)
            
        self.SourceCat=self.SourceCat[self.SourceCat.Exclude==False]

        ClusterList=sorted(list(set(self.SourceCat.Cluster.tolist())))
        # self.NDir=len(ClusterList)
        
        for iCluster,iNewCluster in zip(ClusterList,list(range(self.NDir))):
            ind=np.where(self.SourceCat.Cluster==iCluster)[0]
            self.SourceCat.Cluster[ind]=iNewCluster
            self.REGName=False

        self.Rename()

        self.REGFile=None
        self.MakeREG()


        self.Dirs=sorted(list(set(self.SourceCat.Cluster.tolist())))

        self.WeightDirKeep=np.zeros((self.NDir,),float)
        for diri in self.Dirs:
            ind=np.where(self.SourceCat.Cluster==diri)[0]
            self.WeightDirKeep[diri]=np.sqrt(np.sum(self.SourceCat.Sref[ind]))
        self.WeightDir=self.WeightDirKeep.copy()

        self.ExistToSub=False
        self.ExistToSub=(np.count_nonzero(self.SourceCat.kill==-1)>0)

        self.BuildClusterCat()


    def Rename(self):
        for diri in range(self.NDir):
            ind=np.where(self.SourceCat.Cluster==diri)[0]
            #CatSel=self.SourceCat[self.SourceCat.Cluster==diri]
            Names=["c%is%i."%(diri,i) for i in range(ind.shape[0])]
            self.SourceCat.Name[ind]=Names
        self.REGName=True



    def cluster(self,nk=10,DoPlot=False,PreClusterCat=None,FromClusterCat=""):

        # pylab.clf()
    
        #s.fill(0.)
        #s[0]=1

        cos=np.cos
        sin=np.sin


        self.SourceCat.Cluster=-1
        indSubSel=np.arange(self.SourceCat.shape[0])
        NPreCluster=0

        # #######################################
        # if (PreClusterCat is not None)&(FromClusterCat==""):
        #     N=PreClusterCat.shape[0]
        #     Ns=self.SourceCat.ra.shape[0]
        #     for iReg in range(N):
        #         #d=np.sqrt((self.SourceCat.ra-PreClusterCat.ra[iReg])**2+(self.SourceCat.dec-PreClusterCat.dec[iReg])**2)
        #         ra1=self.SourceCat.ra
        #         ra2=PreClusterCat.ra[iReg]
        #         d1=self.SourceCat.dec
        #         d2=PreClusterCat.dec[iReg]
        #         cosD = sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(ra1-ra2)
        #         d=np.arccos(cosD)
        #         self.SourceCat.Cluster[d<PreClusterCat.Radius[iReg]]=PreClusterCat.Cluster[iReg]
        #         self.SourceCat.Exclude[d<PreClusterCat.Radius[iReg]]=False
        #     print self.SourceCat.Cluster
        #     indPreCluster=np.where(self.SourceCat.Cluster!=-1)[0]
        #     NPreCluster=np.max(PreClusterCat.Cluster)+1
        #     SourceCatPreCluster=self.SourceCat[indPreCluster]
        #     indSubSel=np.where(self.SourceCat.Cluster==-1)[0]
        #     print "number of preselected clusters: %i"%NPreCluster
        
        # if nk==-1:
        #     print "Removing non-clustered sources"%NPreCluster
        #     self.SourceCat=self.SourceCat[self.SourceCat.Cluster!=-1]
        #     print self.SourceCat.Cluster
        #     return
        # #######################################

        SourceCat=self.SourceCat[indSubSel]

        indIMax=np.argmax(SourceCat.I)

        # This is already defined in __init__
        if self.rarad is None:
            self.rarad=SourceCat.ra[indIMax]#np.sum(SourceCat.I*SourceCat.ra)/np.sum(SourceCat.I)
            self.decrad=np.sum(SourceCat.I*SourceCat.dec)/np.sum(SourceCat.I)
            self.CoordMachine = ModCoord.ClassCoordConv(self.rarad, self.decrad)
        
        x,y,s=SourceCat.ra,SourceCat.dec,SourceCat.I
        x,y=self.radec2lm_scalar(x,y)

        # import pylab
        # pylab.clf()
        # pylab.scatter(x,y)
        # pylab.draw()
        # pylab.show()

        
        SourceCat.Cluster=0

        if self.ClusterMethod==1:
            CM=ClassClusterClean(x,y,s,nk,DoPlot=DoPlot)
        elif self.ClusterMethod==2:
            CM=ClassClusterTessel(x,y,s,nk,DoPlot=DoPlot)
        elif self.ClusterMethod==3:
            CM=ClassClusterRadial(x,y,s,nk,DoPlot=DoPlot)
        elif self.ClusterMethod==4:
            if PreClusterCat is not None:
                l0,m0=self.radec2lm_scalar(PreClusterCat.ra,PreClusterCat.dec)
                CM=ClassClusterKMean(x,y,s,nk,DoPlot=DoPlot,PreCluster=(l0,m0))
            else:
                CM=ClassClusterClean(x,y,s,nk,DoPlot=0)#DoPlot)
                DictNode=CM.Cluster()
                ra,dec=[],[]
                ListL,ListM=[],[]
                for idDir in list(DictNode.keys()):
                    LDir=DictNode[idDir]["ListCluster"]
                    #ra0,dec0=np.mean(self.SourceCat.ra[LDir]),np.mean(self.SourceCat.dec[LDir])
                    # print idDir,ra0,dec0,self.SourceCat.ra[LDir].min(),self.SourceCat.ra[LDir].max()
                    # if not np.isnan(ra0):
                    #     ra.append(ra0)
                    #     dec.append(dec0)
                    This_l,This_m=self.radec2lm_scalar(np.array(self.SourceCat.ra[LDir]),np.array(self.SourceCat.dec[LDir]))
                    ListL.append(np.mean(This_l))
                    ListM.append(np.mean(This_m))

                # l0,m0=self.radec2lm_scalar(np.array(ra),np.array(dec))
                l0,m0=np.array(ListL),np.array(ListM)
                nk=l0.size

                CM=ClassClusterKMean(x,y,s,nk,DoPlot=DoPlot,InitLM=(l0,m0))
        elif self.ClusterMethod==5:
            CM=ClassClusterSquareRadial(x,y,s,nk,DoPlot=DoPlot,D_FITS=self.D_FITS)


        REGFile="%s.tessel.reg"%self.TargetList

        #print FromClusterCat
        if FromClusterCat=="":
            DictNode=CM.Cluster()
        else:

            DictNode={}
            if not("reg" in FromClusterCat):
                SourceCatRef=np.load(FromClusterCat)
                SourceCatRef=SourceCatRef.view(np.recarray)
            else:
                R=ModRegFile.RegToNp(FromClusterCat)
                R.Read()
                SourceCatRef=R.CatSel

            ClusterList=sorted(list(set(SourceCatRef.Cluster.tolist())))
            xc,yc=self.radec2lm_scalar(SourceCatRef.ra,SourceCatRef.dec)
            lc=np.zeros((len(ClusterList),),dtype=np.float32)
            mc=np.zeros((len(ClusterList),),dtype=np.float32)
            # for iCluster in ClusterList:
            #     indC=np.where(SourceCatRef.Cluster==iCluster)[0]
            #     lc[iCluster]=np.sum(SourceCatRef.I[indC]*xc[indC])/np.sum(SourceCatRef.I[indC])
            #     mc[iCluster]=np.sum(SourceCatRef.I[indC]*yc[indC])/np.sum(SourceCatRef.I[indC])
            lc,mc=xc,yc
            Ns=x.size
            
            Nc=lc.size
            D=np.sqrt((x.reshape((Ns,1))-lc.reshape((1,Nc)))**2+(y.reshape((Ns,1))-mc.reshape((1,Nc)))**2)
            Cid=np.argmin(D,axis=1)

            #pylab.clf()
            for iCluster in ClusterList:
                ind=np.where(Cid==iCluster)[0]
                DictNode["%3.3i"%iCluster]={}
                DictNode["%3.3i"%iCluster]["ListCluster"]=ind.tolist()
                # pylab.scatter(x[ind],y[ind],c=np.ones((ind.size,))*iCluster,vmin=0,vmax=Nc,lw=0)
                # pylab.draw()
                # pylab.show(block=False)
            


        try:
            CM.ToReg(REGFile,self.rarad,self.decrad)
        except:
            pass

        iK=NPreCluster
        self.NDir=len(list(DictNode.keys()))
        
        # print self.SourceCat.Cluster.min(),self.SourceCat.Cluster.max()
        
        for key in list(DictNode.keys()):
            ind=np.array(DictNode[key]["ListCluster"])
            if ind.size==0: 
                log.print("Direction %i is empty"%int(key))
                iK+=1
                continue
            self.SourceCat["Cluster"][indSubSel[ind]]=iK
            iK+=1



        # if PreClusterCat is not None:
        #     SourceCat=np.concatenate((SourceCatPreCluster,SourceCat))
        #     SourceCat=SourceCat.view(np.recarray)

        #print self.SourceCat.Cluster.min(),self.SourceCat.Cluster.max()
        #self.SourceCat=SourceCat

        
    def ComputeClusterCatWeightedPos(self):
        ClusterCatMean=np.zeros((self.ClusterCat.ra.size,),dtype=dtypeClusterCat)
        ClusterCatMean=ClusterCatMean.view(np.recarray)
        
        ClusterCat=self.ClusterCat
        ClusterCatMean.Name[:]=ClusterCat.Name[:]
        ClusterCatMean.SumI[:]=ClusterCat.SumI[:]
        ClusterCatMean.Cluster[:]=ClusterCat.Cluster[:]
        
        SourceCat=self.SourceCat
        Dirs=list(np.unique(ClusterCatMean.Cluster))
        for iDir in Dirs:
            indDirComps=np.where(SourceCat.Cluster==iDir)[0]
            ra=SourceCat.ra[indDirComps]
            dec=SourceCat.dec[indDirComps]
            S=SourceCat.I[indDirComps]
            ram=np.sum(ra*S)/np.sum(S)
            decm=np.sum(dec*S)/np.sum(S)

            indDirComps=np.where(ClusterCatMean.Cluster==iDir)[0]
            #if indDirComps.size==0: continue
            ClusterCatMean.ra[indDirComps[0]]=ram
            ClusterCatMean.dec[indDirComps[0]]=decm
        self.ClusterCatMeanSource=ClusterCatMean

    
        
    def AppendRefSource(self,rac,decc):
        S0=1e-10
        CatCopy=self.SourceCat[0:1].copy()
        CatCopy['Name']="Reference"
        CatCopy['ra']=rac
        CatCopy['dec']=decc
        CatCopy['Sref']=S0
        CatCopy['I']=S0
        CatCopy['Q']=S0
        CatCopy['U']=S0
        CatCopy['V']=S0
        CatCopy['Cluster']=0
        CatCopy['Type']=0
        self.SourceCat.Cluster+=1
        self.SourceCat=np.concatenate([CatCopy,self.SourceCat])
        self.SourceCat=self.SourceCat.view(np.recarray)
        self.NDir+=1
        self.NSources+=1

        CatCopy=self.ClusterCat[0:1].copy()
        CatCopy['Name']="Reference"
        CatCopy['ra']=rac
        CatCopy['dec']=decc
        CatCopy['SumI']=0
        CatCopy['Cluster']=0
        self.ClusterCat.Cluster+=1
        self.ClusterCat=np.concatenate([CatCopy,self.ClusterCat])
        self.ClusterCat=self.ClusterCat.view(np.recarray)

    def AppendFromSMClass(self,SM):
        #self.SourceCat=np.concatenate([self.SourceCat, SM.SourceCat])
        #self.SourceCat=self.SourceCat.view(np.recarray)
        ClusterCatCopy=SM.ClusterCat.copy()
        ClusterCatCopy.Cluster=ClusterCatCopy.Cluster+np.max(self.ClusterCat.Cluster)+1
        self.ClusterCat=np.concatenate([self.ClusterCat, ClusterCatCopy])
        self.ClusterCat=self.ClusterCat.view(np.recarray)
        
        self.NDir=self.ClusterCat.shape[0]
        self.NSources=self.SourceCat.shape[0]

    def SaveNP(self):
        infile=self.NpFile
        np.save(self.NpFile,self.SourceCat)
        print("   - Numpy catalog file: %s"%ModColor.Str(self.NpFile,col="green"))

    def SelectSourceMouse(self):
        from ClassSelectMouse2 import ClassSelectMouse
        M=ClassSelectMouse()
        ra=self.SourceCat.ra*180/np.pi
        dec=self.SourceCat.dec*180/np.pi
        ra-=np.mean(ra)
        dec-=np.mean(dec)
        M.DefineXY((ra,dec),np.log10(self.SourceCat.I))
        self.SourceCat.Select=M.Start()



    def CutEmptyDirs(self):

        self.ClusterCatOrig=self.ClusterCat.copy()
        # there can be missing Cluster number in SourceCat.Cluster
        self.NDirsOrig=self.ClusterCat.shape[0]
        self.NDir=self.ClusterCat.shape[0]

        AppFlux=self.ClusterCat.SumI
        D={}
        iDirNew=0
        Keep=np.zeros((self.NDir,),bool)
        HasRemoved=0
        
        for iDir in range(self.ClusterCat.size):
            if AppFlux[iDir]==0:
                HasRemoved=1
            else:
                D[iDirNew]=iDir
                iDirNew+=1
                Keep[iDir]=1

        self.MapClusterCatOrigToCut=Keep
        self.ClusterCat=self.ClusterCat[Keep].copy()
        NonZeroDirs=(np.where(Keep)[0]).tolist()
        
        for iClusterNew,iClusterOrig in enumerate(NonZeroDirs):
            ind=np.where(self.ClusterCat.Cluster==iClusterOrig)[0]
            self.ClusterCat.Cluster[ind]=iClusterNew
            ind=np.where(self.SourceCat.Cluster==iClusterOrig)[0]
            self.SourceCat.Cluster[ind]=iClusterNew
            
        self.Dirs=np.sort(np.unique(self.SourceCat.Cluster))#arange(self.ClusterCat.size).tolist()#sorted((np.where(Keep==1)[0]).tolist())
        #self.Dirs=sorted((np.where(Keep==1)[0]).tolist())
        self.NDir=len(self.Dirs)
        if not HasRemoved:
            log.print(ModColor.Str("All directions have been kept in the solve"))
        else:
            log.print(ModColor.Str("There are %i directions left in the solve"%self.NDir))

        
            
        
        
    def BuildClusterCat(self):
        if self.ClusterCatExt is not None:
            self.ClusterCat=self.ClusterCatExt.copy()
            self.ClusterCatOrig=self.ClusterCatExt.copy()
            return
        ClusterCat=np.zeros((len(self.Dirs),),dtype=dtypeClusterCat)
        ClusterCat=ClusterCat.view(np.recarray)
        icat=0
        for d in self.Dirs:
            cat=self.SourceCat[self.SourceCat.Cluster==d]
            l,m=self.CoordMachine.radec2lm(cat.ra, cat.dec)
            if cat.I.max()>0:
                lmean,mmean=np.sum(l*cat.I)/np.sum(cat.I),np.sum(m*cat.I)/np.sum(cat.I)
            else:
                lmean,mmean=np.mean(l),np.mean(m)
            ramean,decmean=self.CoordMachine.lm2radec(np.array([lmean]),np.array([mmean]))
            ClusterCat.ra[icat]=ramean
            ClusterCat.dec[icat]=decmean
            ClusterCat.SumI[icat]=np.sum(cat.I)
            ClusterCat.Cluster[icat]=d
            if "Patch" in cat.dtype.fields.keys() and len(cat["Patch"][0].decode("ascii"))>0:
                ClusterCat.Name[icat]=cat["Patch"][0].decode("ascii")
            else:
                ClusterCat.Name[icat]=cat["Name"][0].decode("ascii").split(".")[-1]

            icat+=1

        #print ClusterCat.ra
        self.ClusterCat=ClusterCat
        self.ClusterCatOrig=self.ClusterCat.copy()


    def radec2lm_scalar(self,ra,dec,rarad0=None,decrad0=None):
        if rarad0==None:
            rarad0=self.rarad
            decrad0=self.decrad
        l = np.cos(dec) * np.sin(ra - rarad0)
        m = np.sin(dec) * np.cos(decrad0) - np.cos(dec) * np.sin(decrad0) * np.cos(ra - rarad0)
        return l,m





    def Calc_LM(self,rac,decc):
        Cat=self.SourceCat
        if not("l" in list(Cat.dtype.fields.keys())):
            Cat=RecArrayOps.AppendField(Cat,'l',float)
            Cat=RecArrayOps.AppendField(Cat,'m',float)
        Cat.l[:],Cat.m[:]=self.radec2lm_scalar(self.SourceCat.ra,self.SourceCat.dec,rac,decc)
        self.SourceCat=Cat
        self.SourceCatKeepForSelector=self.SourceCat.copy()

        Cat=self.ClusterCat
        if not("l" in list(Cat.dtype.fields.keys())):
            Cat=RecArrayOps.AppendField(Cat,'l',float)
            Cat=RecArrayOps.AppendField(Cat,'m',float)
        Cat.l,Cat.m=self.radec2lm_scalar(self.ClusterCat.ra,self.ClusterCat.dec,rac,decc)
        self.ClusterCat=Cat

        Cat=self.ClusterCatOrig
        if not("l" in list(Cat.dtype.fields.keys())):
            Cat=RecArrayOps.AppendField(Cat,'l',float)
            Cat=RecArrayOps.AppendField(Cat,'m',float)
        Cat.l,Cat.m=self.radec2lm_scalar(Cat.ra,Cat.dec,rac,decc)
        self.ClusterCatOrig=Cat



        
    # def Calc_LM(self,rac,decc):
    #     Cat=self.SourceCat
    #     if not("l" in Cat.dtype.fields.keys()):
    #         Cat=RecArrayOps.AppendField(Cat,('l',float))
    #         Cat=RecArrayOps.AppendField(Cat,('m',float))
    #     Cat.l,Cat.m=self.radec2lm_scalar(self.SourceCat.ra,self.SourceCat.dec,rac,decc)
    #     self.SourceCat=Cat
    #     self.SourceCatKeepForSelector=self.SourceCat.copy()
        

    def MakeREG(self):
        self.REGFile="%s.reg"%self.TargetList
        f=open(self.REGFile,"w")
        self.REGName=True
        f.write("# Region file format: DS9 version 4.1\n")
        f.write('global color=green dashlist=8 3 width=1 font="helvetica 7 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        for i in range(self.SourceCat.shape[0]):
            #ss="fk5;ellipse(213.202544,49.871826,0.003909,0.003445,181.376961) # text={P1C1}"
            ra=self.SourceCat.ra[i]*180./np.pi
            dec=self.SourceCat.dec[i]*180./np.pi
            Type=self.SourceCat.Type[i]
            Gmaj=self.SourceCat.Gmaj[i]*180./np.pi*(2.*np.sqrt(2.*np.log(2)))
            Gmin=self.SourceCat.Gmin[i]*180./np.pi*(2.*np.sqrt(2.*np.log(2)))
            if Gmin==0.: Gmin=1./3600
            PA=(self.SourceCat.Gangle[i]+np.pi/2.)*180./np.pi
            rad=20./2600

            #ss="fk5;ellipse(%f,%f,%f,%f,%f) # text={%s}"%(ra,dec,Gmaj,Gmin,0,str(i))
            if self.REGName:
                if Type==1:
                    ss="fk5;ellipse(%f,%f,%f,%f,%f) # text={%s} color=green width=2 "%(ra,dec,Gmaj,Gmin,PA,self.SourceCat.Name[i])
                else:
                    ss="fk5;point(%f,%f) # text={%s} point=circle 5 color=red width=2"%(ra,dec,self.SourceCat.Name[i])
            else:
                if Type==1:
                    ss="fk5;ellipse(%f,%f,%f,%f,%f) # color=green width=2 "%(ra,dec,Gmaj,Gmin,PA)
                else:
                    ss="fk5;point(%f,%f) # point=circle 5 color=red width=2"%(ra,dec)
                

            f.write(ss+"\n")
        f.close()

    def RestoreCat(self):
        self.SourceCat=self.SourceCatKeepForSelector.copy()
        self.Dirs=sorted(list(set(self.SourceCat.Cluster.tolist())))
        self.NDir=len(self.Dirs)
        self.NSources=self.SourceCat.shape[0]
        selDir=np.array(sorted(list(set(self.SourceCat.Cluster.tolist()))))
        # self.WeightDir=self.WeightDirKeep[selDir].copy()


    def SelectSubCat(self,Selector):
        self.Selector=Selector
        self.SourceCat=(self.SourceCatKeepForSelector[self.Selector]).copy()
        self.Dirs=sorted(list(set(self.SourceCat.Cluster.tolist())))
        self.NDir=len(self.Dirs)
        self.NSources=self.SourceCat.shape[0]
        selDir=np.array(sorted(list(set(self.SourceCat.Cluster.tolist()))))
        # self.WeightDir=self.WeightDirKeep[selDir].copy()


        





    def print_sm2(self):
        CatIn=self.SourceCat
        ind=np.argsort(CatIn.Cluster)
        Cat=CatIn[ind]
        TEMPLATE = ('  %(Cluster)5s %(name)10s %(RA)15s %(DEC)15s %(Flux)10s %(alpha)10s %(RefFreq)10s %(Kill)6s ')
        print()
        
        print(" TARGET LIST: ")
        print(TEMPLATE % {
                'Cluster': "K".center(5),
                'name': "Name".center(10),
                'RA': "RA".center(15),
                'DEC': "DEC".center(15),
                'Flux': "Flux".rjust(10),
                'alpha': "alpha".rjust(10),
                'RefFreq': "RefFreq".rjust(10),
                'Kill': "Kill" })

        for i in range(Cat.shape[0]):
            SName=Cat.Name[i]
            SRa=rad2hmsdms.rad2hmsdms(Cat.ra[i],Type="ra").replace(" ",":")
            SDec=rad2hmsdms.rad2hmsdms(Cat.dec[i]).replace(" ",".")
            SI="%6.3f"%Cat.I[i]
            SAlpha="%4.2f"%Cat.alpha[i]
            SRefFreq="%5.1f"%(Cat.RefFreq[i]/1.e6)
            SKill="%i"%Cat.kill[i]
            SCluster="%2.2i"%Cat.Cluster[i]
            StrOut = TEMPLATE % {
                'Cluster': SCluster.center(5),
                'name': SName.center(10),
                'RA': SRa,
                'DEC': SDec,
                'Flux': SI,
                'alpha': SAlpha,
                'RefFreq':SRefFreq,
                'Kill':SKill }
            print(StrOut)
                

    def print_sm(self,Cat):
        if self.infile_cluster=="":
            print(" TARGET LIST: ")
            format="%13s%20s%20s%10s%10s%10s"#%10s"
            print(format%("Name","Ra","Dec","Flux","alpha","RefFreq"))#,"Kill")
            for i in range(Cat.shape[0]):
    
    
                SName=Cat.Name[i]
                SRa=rad2hmsdms.rad2hmsdms(Cat.ra[i]/15).replace(" ",":")
                SDec=rad2hmsdms.rad2hmsdms(Cat.dec[i]).replace(" ",".")
                SI=Cat.I[i]
                SAlpha=Cat.alpha[i]
                SRefFreq=Cat.RefFreq[i]
                SKill=str(Cat.kill[i]==1)
                #print "%13s%20s%20s%10.4f%10s%10.2e%8s"%(SName,SRa,SDec,SI,SAlpha,SRefFreq,SKill)
                print("%13s%20s%20s%10.4f%10s%10.2e"%(SName,SRa,SDec,SI,SAlpha,SRefFreq))#,SKill)
        else:
            format="%10s%10s%15s%15s%10s%10s%10s"#%8s"
            print()
            print(" TARGET LIST: ")
            print(format%("Group","Name","Ra","Dec","Flux","alpha","RefFreq"))#,"kill")
            for i in range(np.max(Cat.Cluster)+1):
                ind=np.where(Cat.Cluster==i)[0]
                for j in range(ind.shape[0]):
                    jj=ind[j]
                    SName=Cat.Name[jj]
                    SRa=rad2hmsdms.rad2hmsdms(Cat.ra[jj]/15).replace(" ",":")
                    SDec=rad2hmsdms.rad2hmsdms(Cat.dec[jj]).replace(" ",".")
                    SI=Cat.I[jj]
                    SAlpha=Cat.alpha[jj]
                    SRefFreq=Cat.RefFreq[jj]
                    SKill=str(Cat.kill[jj]==1)
                    SGroup=str(i)
                    #print "%10s%10s%15s%15s%8.4f%10s%10.2e%8s"%(SGroup,SName,SRa,SDec,SI,SAlpha,SRefFreq)#,SKill)
                    print("%10s%10s%15s%15s%8.4f%10s%10.2e"%(SGroup,SName,SRa,SDec,SI,SAlpha,SRefFreq))#,SKill)
                    

    
    
    
    
