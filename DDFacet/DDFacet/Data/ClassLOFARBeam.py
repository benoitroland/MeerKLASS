'''
DDFacet, a facet-based radio imaging package
Copyright (C) 2013-2016  Cyril Tasse, l'Observatoire de Paris,
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
'''

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from DDFacet.compatibility import range

from DDFacet.Other import logger
log= logger.getLogger("ClassLOFARBeam")
from DDFacet.Other import ClassTimeIt
from DDFacet.Other import ModColor
#import everybeam
from astropy.coordinates import AltAz, EarthLocation, ITRS, SkyCoord
from astropy.time import Time
import astropy.units as u

import numpy as np

def test():
    from DDFacet.Data import ClassMS
    GD={"Beam":{"PhasedArrayMode":"AE",
                "DtBeamMin":5,
                "NBand":2},
        "Cache":{"Dir":"."}}
    MS=ClassMS.ClassMS("L526161_SB347_uv.pre-cal_1288A51A4t_142MHz.pre-cal.ms.send.tsel",DoReadData=0)
    rac,decc=MS.radec
    dd=3*np.pi/180
    N=11
    ra,dec=np.mgrid[rac-dd:rac+dd:N*1j,decc-dd:decc+dd:N*1j]
    ra=ra.ravel()
    dec=dec.ravel()
    time=MS.F_tstart+4*3600

    from DDFacet.Other import ClassTimeIt

    T=ClassTimeIt.ClassTimeIt()
    BeamMachine=ClassLOFARBeam(MS,GD)
    J0=BeamMachine.GiveInstrumentBeam(time,ra,dec)
    T.timeit("J0")
    
    BeamMachine=ClassLOFARBeamEB(MS,GD)
    J1=BeamMachine.GiveInstrumentBeam(time,ra,dec)
    T.timeit("J1")
    
    nd,na,nf,_,_=J0.shape
    J0=J0.reshape((N,N,na,nf,2,2))
    J1=J1.reshape((N,N,na,nf,2,2))

    iAnt=10
    op0=lambda x: np.log10(np.abs(x))
    op0=lambda x: np.abs(x)
    op1=np.angle
    import pylab
    pylab.clf()
    pylab.subplot(2,2,1)
    pylab.imshow(op0(J0[:,:,iAnt,0,0,0]),vmin=0,vmax=1,interpolation="nearest")
    pylab.subplot(2,2,2)
    pylab.imshow(op0(J1[:,:,iAnt,0,0,0]),vmin=0,vmax=1,interpolation="nearest")
    pylab.subplot(2,2,3)
    pylab.imshow(op1(J0[:,:,iAnt,0,0,0]),vmin=-np.pi,vmax=np.pi,interpolation="nearest")
    pylab.subplot(2,2,4)
    pylab.imshow(op1(J1[:,:,iAnt,0,0,0]),vmin=-np.pi,vmax=np.pi,interpolation="nearest")
    pylab.draw()
    pylab.show()




class ClassLOFARBeamEB():
    def __init__(self,MS,GD):
        self.GD=GD
        self.MS=MS
        self.SR=None
        self.InitLOFARBeam()
        self.CalcFreqDomains()

    def InitLOFARBeam(self):
        GD=self.GD
        PhasedArrayMode=GD["Beam"]["PhasedArrayMode"]
        print("  LOFAR beam model in %s mode"%(PhasedArrayMode), file=log)
        self.useArrayFactor=("A" in PhasedArrayMode)
        self.useElementBeam=("E" in PhasedArrayMode)
        if self.SR is not None: return
        import everybeam
        use_differential_beam=1#(self.useArrayFactor and not self.useElementBeam)
        self.SR = everybeam.load_telescope(self.MS.MSName)#,
        #use_differential_beam=use_differential_beam,
        #                                   use_channel_frequency=1)#use_differential_beam)
        self.SReval=self.SR.station_response
        
        
    def getBeamSampleTimes(self,times, **kwargs):
        DtBeamMin = self.GD["Beam"]["DtBeamMin"]
        DtBeamSec = DtBeamMin*60
        tmin=times[0]
        tmax=times[-1]+1
        TimesBeam=np.arange(tmin,tmax,DtBeamSec).tolist()
        if not(tmax in TimesBeam): TimesBeam.append(tmax)
        return TimesBeam

    def getFreqDomains(self):
        return self.FreqDomains

    def CalcFreqDomains(self):
        ChanWidth=self.MS.ChanWidth.ravel()[0]
        ChanFreqs=self.MS.ChanFreq.flatten()

        NChanJones=self.GD["Beam"]["NBand"]
        if NChanJones==0:
            NChanJones=self.MS.NSPWChan
        ChanEdges=np.linspace(ChanFreqs.min()-ChanWidth/2.,ChanFreqs.max()+ChanWidth/2.,NChanJones+1)

        FreqDomains=[[ChanEdges[iF],ChanEdges[iF+1]] for iF in range(NChanJones)]
        FreqDomains=np.array(FreqDomains)
        self.FreqDomains=FreqDomains
        self.NChanJones=NChanJones

        MeanFreqJonesChan=(FreqDomains[:,0]+FreqDomains[:,1])/2.
        DFreq=np.abs(self.MS.ChanFreq.reshape((self.MS.NSPWChan,1))-MeanFreqJonesChan.reshape((1,NChanJones)))
        self.VisToJonesChanMapping=np.argmin(DFreq,axis=1)

    def GiveRawBeam(self,time,ra,dec):
        #self.LoadSR()
        Beam=np.zeros((ra.shape[0],self.MS.na,self.MS.NSPWChan,2,2),dtype=np.complex)
        for i in range(ra.shape[0]):
            for iAnt in range(self.MS.na):
                for freq in self.MS.ChanFreq.flatten():
                    Beam[i]=self.SReval(time,iAnt,freq,ra[i],dec[i])
                    #print(time,iAnt,freq,ra[i],dec[i])
        #Beam=np.swapaxes(Beam,1,2)
        return Beam

    def GiveInstrumentBeam(self,*args,**kwargs):
        
        T=ClassTimeIt.ClassTimeIt("GiveInstrumentBeam")
        T.disable()
        Beam=self.GiveRawBeam(*args,**kwargs)
        nd,na,nch,_,_=Beam.shape
        T.timeit("0")
        MeanBeam=np.zeros((nd,na,self.NChanJones,2,2),dtype=Beam.dtype)
        for ich in range(self.NChanJones):
            indCh=np.where(self.VisToJonesChanMapping==ich)[0]
            MeanBeam[:,:,ich,:,:]=np.mean(Beam[:,:,indCh,:,:],axis=2)
        T.timeit("1")

        return MeanBeam


class ClassLOFARBeam():
    def __init__(self,MS,GD):
        self.GD=GD
        self.MS=MS
        self.SR=None
        self.InitLOFARBeam()
        self.CalcFreqDomains()

    def InitLOFARBeam(self):
        GD=self.GD
        PhasedArrayMode=GD["Beam"]["PhasedArrayMode"]
        print("  LOFAR beam model in %s mode"%(PhasedArrayMode), file=log)
        useArrayFactor=("A" in PhasedArrayMode)
        useElementBeam=("E" in PhasedArrayMode)
        if self.SR is not None: return

        import lofar.stationresponse as lsr

        self.SR = lsr.stationresponse(self.MS.MSName,
                                      useElementResponse=useElementBeam,
                                      #useElementBeam=useElementBeam,
                                      useArrayFactor=useArrayFactor,
                                      useChanFreq=True)
        self.SR.setDirection(self.MS.rarad,self.MS.decrad)


    def getBeamSampleTimes(self,times, **kwargs):
        DtBeamMin = self.GD["Beam"]["DtBeamMin"]
        DtBeamSec = DtBeamMin*60
        tmin=times[0]
        tmax=times[-1]+1
        TimesBeam=np.arange(tmin,tmax,DtBeamSec).tolist()
        if not(tmax in TimesBeam): TimesBeam.append(tmax)
        return TimesBeam

    def getFreqDomains(self):
        return self.FreqDomains

    def CalcFreqDomains(self):
        ChanWidth=self.MS.ChanWidth.ravel()[0]
        ChanFreqs=self.MS.ChanFreq.flatten()

        NChanJones=self.GD["Beam"]["NBand"]
        if NChanJones==0:
            NChanJones=self.MS.NSPWChan
        ChanEdges=np.linspace(ChanFreqs.min()-ChanWidth/2.,ChanFreqs.max()+ChanWidth/2.,NChanJones+1)

        FreqDomains=[[ChanEdges[iF],ChanEdges[iF+1]] for iF in range(NChanJones)]
        FreqDomains=np.array(FreqDomains)
        self.FreqDomains=FreqDomains
        self.NChanJones=NChanJones

        MeanFreqJonesChan=(FreqDomains[:,0]+FreqDomains[:,1])/2.
        DFreq=np.abs(self.MS.ChanFreq.reshape((self.MS.NSPWChan,1))-MeanFreqJonesChan.reshape((1,NChanJones)))
        self.VisToJonesChanMapping=np.argmin(DFreq,axis=1)

    def GiveRawBeam(self,time,ra,dec):
        #self.LoadSR()
        Beam=np.zeros((ra.shape[0],self.MS.na,self.MS.NSPWChan,2,2),dtype=np.complex)
        for i in range(ra.shape[0]):
            self.SR.setDirection(ra[i],dec[i])
            Beam[i]=self.SR.evaluate(time)
        #Beam=np.swapaxes(Beam,1,2)
        return Beam

    def GiveInstrumentBeam(self,*args,**kwargs):
        
        T=ClassTimeIt.ClassTimeIt("GiveInstrumentBeam")
        T.disable()
        Beam=self.GiveRawBeam(*args,**kwargs)
        nd,na,nch,_,_=Beam.shape
        T.timeit("0")
        MeanBeam=np.zeros((nd,na,self.NChanJones,2,2),dtype=Beam.dtype)
        for ich in range(self.NChanJones):
            indCh=np.where(self.VisToJonesChanMapping==ich)[0]
            MeanBeam[:,:,ich,:,:]=np.mean(Beam[:,:,indCh,:,:],axis=2)
        T.timeit("1")

        return MeanBeam


