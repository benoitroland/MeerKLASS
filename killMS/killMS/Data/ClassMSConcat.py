import numpy as np
from .ClassMS import ClassMS
from DDFacet.Other import logger
log = logger.getLogger("ClassMSConcat")
from killMS.Other import ModColor
from killMS.Other.rad2hmsdms import rad2hmsdms
from killMS.Other import ModColor

class ClassMSConcat():
    def __init__(self,*args,**kwargs):
        self.kwargs=kwargs
        self.args=args
        self.Type="Virtual"

        self.MSName=MSName=args[0]
        MSName0 = MSName
        ListMSName = [ l.strip() for l in open(MSName).readlines() ]
        self.MSList=ListMSName
        self.ListMS=[ClassMS(MSName,**kwargs) for MSName in self.MSList]
        
        self.MS0=self.ListMS[0]
        for key in self.MS0.__dict__.keys():
            if key not in ["MSName"]:
                setattr(self,key,self.MS0.__dict__[key])

        CheckKey=['nbl', 'F_nrows', 'F_A0', 'F_A1', 'na', 'F_tstart', 'F_times_all', 'F_times', 'F_ntimes', 'dt', 'tMean', 'DTs', 'DTh', 'radec', 'OriginalRadec', 'rarad', 'decrad', 'StationNames', 'rac', 'decc']
        log.print("Checking if the %i MSs are aligned for frequency concatenation"%len(self.MSList))

        ThisFreq0=self.MS0.reffreq
        for MS in self.ListMS[1:]:
            ThisFreq=MS.reffreq
            if ThisFreq>ThisFreq0:
                ThisFreq0=ThisFreq
            else:
                raise ValueError("MS must be sorted in increasing freqs")
                
            for key in CheckKey:
                Obj0=getattr(self,key)
                Obj1=getattr(MS,key)
                if type(Obj0).__name__=="ndarray":
                    Test=np.allclose(Obj0,Obj1)
                else:
                    Test=(Obj0==Obj1)
                if not Test:
                    raise ValueError("Invalid MS %s field %s" % (MS.MSName, key))
        log.print("  ... MSs are aligned")


        
        ChanWidth=[]
        ChanFreq=[]
        wavelength_chan=[]
        i0=0
        ListCh0Ch1=[]
        for MS in self.ListMS:
            i1=i0+MS.NSPWChan
            ChanWidth=np.concatenate([ChanWidth,MS.ChanWidth.ravel()])
            ChanFreq=np.concatenate([ChanFreq,MS.ChanFreq.ravel()])
            wavelength_chan=np.concatenate([wavelength_chan,MS.wavelength_chan.ravel()])
            ListCh0Ch1.append((i0,i1))
            i0+=MS.NSPWChan
            
        self.ListCh0Ch1=ListCh0Ch1
        self.NSPWChan=ChanFreq.size
        self.reffreq=np.mean(ChanFreq)
        self.wavelength_chan=wavelength_chan
        self.dFreq=self.ChanWidth=ChanWidth
        
        
        self.ChanFreq=ChanFreq
        
        log.print("There are %i channels in the virtually concatenated MS..."%self.NSPWChan)

    def GiveAllUVW(self,WeightCol="WEIGHT"):
        uvw,W,F=self.MS0.GiveAllUVW(WeightCol=WeightCol)
        for MS in self.ListMS[1:]:
            _,Ws,Fs=MS.GiveAllUVW(WeightCol=WeightCol)
            W=np.concatenate([W,Ws],axis=1)
            F=np.concatenate([F,Fs],axis=1)
        return uvw,W,F

    def LoadSR(self,useElementBeam=True,useArrayFactor=True):
        for MS in self.ListMS:
            MS.LoadSR(useElementBeam=useElementBeam,
                      useArrayFactor=useArrayFactor)

    def GiveBeam(self,time,ra,dec):
        Beam=self.MS0.GiveBeam(time,ra,dec)
        for MS in self.ListMS[1:]:
            Beam = np.concatenate([Beam,MS.GiveBeam(time,ra,dec)],axis=2)
        return Beam

    def AddCol(self,ColName,LikeCol="DATA",ColDesc=None,ColDescDict=None):
        for MS in self.ListMS:
            MS.AddCol(ColName,LikeCol=LikeCol,ColDesc=ColDesc,ColDescDict=ColDescDict)
        
    def putcol(self,ColName,d,row0,nrows):
        for iMS,MS in enumerate(self.ListMS):
            i0,i1=self.ListCh0Ch1[iMS]
            MS.putcol(ColName,d[:,i0:i1,...],row0,nrows)

    def Rotate(self,DATA,RotateType=["uvw","vis"],Sense="ToTarget",DataFieldName="data"):
        for iMS,MS in enumerate(self.ListMS):
            i0,i1=self.ListCh0Ch1[iMS]
            ThisDATA={"uvw":DATA["uvw"],
                      "data":DATA["data"][:,i0:i1,...].copy()}
            MS.Rotate(ThisDATA,
                      RotateType=RotateType,
                      Sense=Sense,
                      DataFieldName=DataFieldName)
            DATA["data"][:,i0:i1,...]=ThisDATA["data"][...]

            
    def ToOrigFreqOrder(self,data):
        
        if self.DoRevertChans:
            dOut=np.zeros_like(data)
            for iMS,MS in enumerate(self.ListMS):
                i0,i1=self.ListCh0Ch1[iMS]
                dOut[:,i0:i1]=data[:,i0:i1][:,::-1]
        else:
            dOut=data
            
        if self.MS0.NPolOrig==2:
            dOut=dOut[:,:,0::3]

        return dOut[:,self.MS0.ChanSlice,...]

    def ReadData(self,t0=0,t1=-1,DoPrint=False,ReadWeight=False):

        HasEmptyChunk=False
        for MS in self.ListMS:
            rep=MS.ReadData(t0=t0,t1=t1,DoPrint=DoPrint,ReadWeight=ReadWeight)
            if rep is None:
                log.print("[%s] This data chunk is empty"%MS.MSName)

        if HasEmptyChunk:
            return None
        
        self.CorrelationNames=self.MS0.CorrelationNames
        self.CurrentChunkTimeRange_SinceT0_sec=self.MS0.CurrentChunkTimeRange_SinceT0_sec
        self.ROW0=self.MS0.ROW0
        self.ROW1=self.MS0.ROW1
        self.nRowRead=self.MS0.nRowRead
        self.Time0=self.MS0.Time0
        self.HasWeights=self.MS0.HasWeights

        self.Weights=self.MS0.Weights
        log.print("Concatenate weight tables...")
        for MS in self.ListMS[1:]:
            self.Weights=np.concatenate([self.Weights,MS.Weights],axis=1)
            
        self.HasWeights=self.MS0.HasWeights
        self.multidata=self.MS0.multidata
        self.ReverseAntOrder=self.MS0.ReverseAntOrder
        self.swapped=self.MS0.swapped
        self.TimeInterVal=self.MS0.TimeInterVal
        self.NPolOrig=self.MS0.NPolOrig

        self.uvw=self.MS0.uvw
        self.uvw_dt=self.MS0.uvw_dt
        self.swapped=self.MS0.swapped
        self.times_all=self.MS0.times_all
        self.times=self.MS0.times
        self.ntimes=self.MS0.ntimes
        self.nrows=self.MS0.nrows
        self.IndFlag=self.MS0.IndFlag
        self.A0=self.MS0.A0
        self.A1=self.MS0.A1

        
        self.data=self.MS0.data
        self.flag_all=self.MS0.flag_all
        log.print("Concatenate data/flag tables...")
        for MS in self.ListMS[1:]:
            self.data=np.concatenate([self.data,MS.data],axis=1)
            self.flag_all=np.concatenate([self.flag_all,MS.flag_all],axis=1)
        return True
    
    def __str__(self):
        ll=[]
        ll.append(ModColor.Str(" VIRTUALLY CONCATENATED MS PROPERTIES: "))
        ll.append("   - File Name: %s"%ModColor.Str(self.MSName,col="green"))
        #ll.append("   - Column Name: %s"%ModColor.Str(str(self.ColName),col="green"))
        ll.append("   - Selection: %s"%( ModColor.Str(str(self.TaQL),col="green")))
        ll.append("   - Pointing center: (ra, dec)=(%s, %s) "%(rad2hmsdms(self.rarad,Type="ra").replace(" ",":")\
                                                               ,rad2hmsdms(self.decrad,Type="dec").replace(" ",".")))
        f0,f1=self.ChanFreq.min()/1e6,self.ChanFreq.max()/1e6
        ll.append("   - Frequency = %.2f MHz [ %.2f -> %.2f MHz]"%(self.reffreq/1e6,f0,f1))
        ll.append("   - Wavelength = %5.2f meters"%(np.mean(self.wavelength_chan)))
        ll.append("   - Time bin = %4.1f seconds"%(self.dt))
        ll.append("   - Total Integration time = %6.2f hours"%self.DTh)
        ll.append("   - Number of antenna  = %i"%self.na)
        ll.append("   - Number of baseline = %i"%self.nbl)
        #ll.append("   - Number of SPW = %i"%self.NSPW)
        ll.append("   - Number of channels = %i"%self.NSPWChan)
        
        s=" ".join(["%.2f"%(x/1e6) for x in self.ChanFreq.flatten()])
        #ll.append("   - Chan freqs = %s"%(ListToStr(s.split(" "),Unit="MHz")))
        
        ss="\n".join(ll)+"\n"
        return ss
