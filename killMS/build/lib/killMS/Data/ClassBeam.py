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
from . import ClassMS
from pyrap.tables import table
from DDFacet.Other import logger
log=logger.getLogger("ClassBeam")
from killMS.Array import ModLinAlg

def GiveNXNYPanels(Ns,ratio=800/500):
    nx=int(round(np.sqrt(Ns/ratio)))
    ny=int(nx*ratio)
    if nx*ny<Ns: ny+=1
    return nx,ny


class ClassBeam():
    def __init__(self,MSName,GD,SM,ColName=None):
        self.GD=GD
        self.SM=SM
        self.MSName=MSName#self.GD["VisData"]["MSName"]
        if ColName is not None:
            self.ColName=ColName
        else:
            self.ColName=self.GD["VisData"]["InCol"]
            
        if self.GD is None or not self.GD["VisData"]["ConcatFreq"]:
            from .ClassMS import ClassMS
        else:
            from .ClassMSConcat import ClassMSConcat as ClassMS
            
        self.MS=ClassMS(self.MSName,Col=self.ColName,DoReadData=False,GD=GD)
        self.DtBeamMin=self.GD["Beam"]["DtBeamMin"]

    def GiveMeanBeam(self,NTimes=None):
        log.print( "Calculate mean beam ... ")
        #t=table(self.MSName,ack=False)
        times=self.MS.F_times_all#t.getcol("TIME")
        #t.close()
        
        DicoBeam=self.GiveBeamMeanAllFreq(times,NTimes=NTimes)
        J=DicoBeam["Jones"]
        AbsMean=np.mean(np.abs(J),axis=0)
        J0=J#[...,0,0]

        JJ=J0*J0.conj()
        AbsMean=np.sqrt(np.mean(np.abs(J0*J0.conj()),axis=0))


        # # #########################
        # import pylab
        # import matplotlib.gridspec as gridspec
        # nt,nd,na,nch,_,_=JJ.shape
        # fig=pylab.figure(0,figsize=(13,8))
        # nx,ny=GiveNXNYPanels(na)
        # gs1 = gridspec.GridSpec(2*nx, ny)
        # gs1.update(wspace=0.05, hspace=0.05, left=0.05, right=0.95, bottom=0.05, top=0.95)
    
        # pylab.clf()
        # op0=np.abs
        # op1=np.angle
        # ADir_0=op0(JJ[:,:,:,0,0,0])
        # ADir_1=op1(JJ[:,:,:,0,0,0])
        # vmin,vmax=0,2
    
        # Mean=np.mean(ADir_0)
        # MAD=np.sqrt(np.median((ADir_0-Mean)**2))
        # vmin,vmax=Mean-10*MAD,Mean+10*MAD
    
        # # if self.Mask is not None:
        # #     M=self.Mask[:,:,:,iDir,0,0]
        # #     ADir_0[M==1]=np.nan
        # #     ADir_1[M==1]=np.nan


            
        # iAnt=0
        # for i in range(nx):
        #     for j in range(ny):
        #         if iAnt>=na:continue
        #         #pylab.title(StationNames[iAnt], fontsize=9)
    
        #         A_0=ADir_0[:,:,iAnt]
        #         A_1=ADir_1[:,:,iAnt]
        #         ax = pylab.subplot(gs1[2*i,j])
        #         #ax2 = ax.twinx()
        #         ax.imshow(A_0.T,vmin=vmin,vmax=vmax,interpolation="nearest",aspect='auto',cmap="gray")
        #         ax.set_xticks([])
        #         ax.set_yticks([])
                
        #         ax = pylab.subplot(gs1[2*i+1,j])
        #         #ax2 = ax.twinx()
        #         #A=Sols.G[:,:,iAnt,iDir,0,0]
        #         F=5.
        #         ax.imshow(A_1.T,vmin=-np.pi/F,vmax=np.pi/F,interpolation="nearest",aspect='auto')
        #         #ax.imshow(A_1.T,interpolation="nearest",aspect='auto')
        #         #print iAnt,np.max(np.abs(A_1))
        #         ax.set_xticks([])
        #         ax.set_yticks([])
        #         iAnt+=1
        # #pylab.suptitle('Channel %i (op0[%5.2f, %5.2f])'%(iChan,vmin,vmax))
        # #pylab.tight_layout(pad=3., w_pad=0.5, h_pad=2.0)
        # pylab.draw()
        # pylab.show()
        # #pylab.pause(0.1)
        # #time.sleep(1)
        

        




        
        return AbsMean
    # def SetLOFARBeam(self,LofarBeam):
    #     self.BeamMode,self.DtBeamMin,self.BeamRAs,self.BeamDECs = LofarBeam
    #     log.print( "Set LOFARBeam in %s Mode"%self.BeamMode)
    #     useArrayFactor=("A" in self.BeamMode)
    #     useElementBeam=("E" in self.BeamMode)
    #     self.MS.LoadSR(useElementBeam=useElementBeam,useArrayFactor=useArrayFactor)
    #     self.ApplyBeam=True
        
        
    def GiveBeamMeanAllFreq(self,times,NTimes=None):
        if self.GD["Beam"]["BeamModel"]=="LOFAR":
            useArrayFactor=("A" in self.GD["Beam"]["PhasedArrayMode"])
            useElementBeam=("E" in self.GD["Beam"]["PhasedArrayMode"])
            self.MS.LoadSR(useElementBeam=useElementBeam,useArrayFactor=useArrayFactor)
        elif self.GD["Beam"]["BeamModel"]=="FITS" or self.GD["Beam"]["BeamModel"]=="ATCA" or self.GD["Beam"]["BeamModel"]=="NENUFAR":
            self.MS.LoadDDFBeam()


        
        tmin,tmax=times[0],times[-1]
        if NTimes is None:
            DtBeamSec=self.DtBeamMin*60
            DtBeamMin=self.DtBeamMin
        else:
            DtBeamSec=(tmax-tmin)/(NTimes+1)
            DtBeamMin=DtBeamSec/60
        log.print( "  Update beam [Dt = %3.1f min] ... "%DtBeamMin)
            
        TimesBeam=np.arange(tmin,tmax,DtBeamSec).tolist()
        if not(tmax in TimesBeam): TimesBeam.append(tmax)
        TimesBeam=np.array(TimesBeam)
        T0s=TimesBeam[:-1]
        T1s=TimesBeam[1:]
        Tm=(T0s+T1s)/2.
        RA,DEC=self.SM.ClusterCat.ra,self.SM.ClusterCat.dec
        NDir=RA.size
        Beam=np.zeros((Tm.size,NDir,self.MS.na,self.MS.NSPWChan,2,2),np.complex64)
        for itime in range(Tm.size):
            ThisTime=Tm[itime]
            Beam[itime]=self.MS.GiveBeam(ThisTime,RA,DEC)
    
        ###### Normalise
        rac,decc=self.MS.OriginalRadec
        if self.GD["Beam"]["CenterNorm"]==1:
            for itime in range(Tm.size):
                ThisTime=Tm[itime]
                Beam0=self.MS.GiveBeam(ThisTime,np.array([rac]),np.array([decc]))
                Beam0inv=ModLinAlg.BatchInverse(Beam0)
                nd,_,_,_,_=Beam[itime].shape
                Ones=np.ones((nd, 1, 1, 1, 1),np.float32)
                Beam0inv=Beam0inv*Ones
                Beam[itime]=ModLinAlg.BatchDot(Beam0inv,Beam[itime])
        ###### 

        nt,nd,na,nch,_,_= Beam.shape
        Beam=np.mean(Beam,axis=3).reshape((nt,nd,na,1,2,2))
        
        
        DicoBeam={}
        DicoBeam["t0"]=T0s
        DicoBeam["t1"]=T1s
        DicoBeam["tm"]=Tm
        DicoBeam["Jones"]=Beam
        return DicoBeam
